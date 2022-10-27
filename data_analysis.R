# this is the main analysis script. it depends on the data_processing script for
# tidy data inputs, and the data_analysis_landcover_tables script for
# pre-generated AICc tables of arthropod~landcover models.

# Load packages ####
library(BAMMtools) # for fast Jenks breaks
library(car) # for Anova() on lmer model objects
library(classInt) # for kmeans clustering
library(corrr) # for correlation plots of landcover vars
library(crayon) # for colored terminal outputs
library(data.table) # for rbindlist() to rbind a list of tables and make id col
library(DiagrammeRsvg) # for plotting SEMs
library(effects) # for effects plots
library(emmeans) # for computing SEM marginal means
library(gridExtra) # create multi-panel plots
library(ggeffects) # for easy effects plots
library(ggfortify) # create PCA plots
library(ggiraphExtra) # more easy effects plots
library(ggpmisc) # make ggplots that show R2 value with stat_poly_eq()
library(ggridges) # for ggridges plots
library(grid) # for grobTree() to make text annotations on violin plots
library(gtools) # for mixedsort() to arrange factor levels in vegdata tibble
library(hardhat) # for get_levels to extract factor levels in a tidy way
library(knitr) # for knitting R markdown docs
library(lme4) # for univariate mixed-effects models
library(lmerTest) # for lmer with p values
library(magrittr) # for assignment and exposition pipes
library(MuMIn) # model selection tools
library(mvabund) # for building multivariate mods of insect density
library(piecewiseSEM) # for structural equation modeling
library(plotly) # interactive plots with plotly()
library(rsvg) # for more SEM plots
library(sjPlot) # create effects plots on lmer objects
library(tidyselect) # for peek()
library(tidytext) # sort ggcols after faceting
library(tidyverse) # R packages for data science
library(varhandle) # easily create dummy vars with to.dummy()
library(vegan) # for diversity indices in vegdata
library(webshot) # to capture tab_model output (or other html output) as png

# Define functions ####
# note: this is not the only place functions are defined

# p value formatting function
# not sure how this works. from stackexchange
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {

  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

# define stupid 'reverse paste' function for mapping paste0() over a list
rpaste0 <- function (x,y) {
  paste0(y,x)
}

# import data ####
## read in tidy csvs ####
# arthropod data by plot, taxa in separate cols
data <- read_csv('tidy_data/data.csv', col_types = 'fffffdddddddddddddddddddd')
# arthropod data by plot, all taxa in one col
data_long <- read_csv('tidy_data/data_long.csv', col_types = 'fffffffd')
# landcover class by field, regular classification
landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffdddddddddddd')
# landcover class by field, target alfalfa fields identified manually
landcoverFixed <-
  read_csv('tidy_data/landcover_fixed.csv', col_types = 'ffffdddddddddddd')
# vegetation survey data by survey plot
vegPlots <- read_csv('tidy_data/vegPlots.csv', col_types = 'fffffff')
# vegetation survey data by site (farm)
vegSites <- read_csv('tidy_data/vegSites.csv', col_types = 'f')

## wrangle ####
# make df of log(total abundance) across treatments
mean_density_long <- data_long %>%
  # divide 'Pre-' values by 3 to account for unequal sampling area
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  # create unique id for each location and sampling time
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  # using the new id col, group treatments together
  # note: 'Full' treatment has already been removed.
  # only Pre-, Control, and Sham remain.
  group_by(id) %>%
  # summarize using the mean
  summarize(Mean_Density = mean(Density)) %>%
  # log+1 transform all values
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  # separate id column to restore original columns
  separate(id, c('Site', 'Field', 'Plot', 'Season', 'Taxa')) %>%
  # create watering method column
  mutate(wateringMethod = case_when(Site %in% c('Fallon', 'Lovelock') ~
                                      'Flooding',
                                    Site == 'Minden' ~
                                      'Sprinklers',
                                    Site == 'Yerington' & Field %in% c(2, 3) ~
                                      'Flooding',
                                    Site == 'Yerington' & Field == 1 ~
                                      'Sprinklers'))

# make wide version of same df
mean_density_wide <- mean_density_long %>%
  # create separate columns for each taxon
  pivot_wider(names_from = Taxa, values_from = Mean_Density)

# create margin data
field_margins <- vegPlots %>% filter(type == 'Margin') %>%
  mutate(total_cover = select(., 12:132) %>% rowSums(na.rm = TRUE)) %>%
  group_by(field_id, season) %>%
  summarize(shan = mean(shan),
            rich = mean(rich),
            totalCover = mean(total_cover),
            .groups = 'keep') %>%
  # make join key
  mutate(id = str_replace(field_id, ' ', '0'), Season = season) %>%
  ungroup() %>%
  select(-field_id, -season)

# summarize mean_density across fields
# don't use mean_density_field here:
# must average within fields BEFORE log transformation for greatest accuracy
mean_density_field <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  # log+1 transform all values
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  separate(id, c('Site', 'Field', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density) %>%
  mutate(id = paste0(Site, 0, Field), .before = Site) %>%
  # create watering method column
  mutate(wateringMethod = case_when(Site %in% c('Fallon', 'Lovelock') ~
                                      'Flooding',
                                    Site == 'Minden' ~
                                      'Sprinklers',
                                    Site == 'Yerington' & Field %in% c(2, 3) ~
                                      'Flooding',
                                    Site == 'Yerington' & Field == 1 ~
                                      'Sprinklers'), .after = Season) %>%
  # join margin data
  left_join(field_margins, by = c('id' = 'id', 'Season' = 'Season')) %>%
  # relocate field margin data cols
  relocate(shan, rich, totalCover, .after = wateringMethod)

## prepare landcover data ####
# landcover data - 4 versions: (7- vs 8-class) * (regular vs fixed)
landcover7 <- landcover %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedyWet = class5 + class8,
         water = class11
  ) %>%
  # scale landcover factors here!
  mutate(across(alfalfa:water, ~ as.vector(scale(.))))


#   mutate_at(21:69, ~ as.vector(scale(.)))

landcover8 <- landcover %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedy = class5,
         wet =  class8,
         water = class11
  ) %>%
  # scale landcover factors here!
  mutate(across(alfalfa:water, ~ as.vector(scale(.))))

landcover7Fixed <- landcoverFixed %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedyWet = class5 + class8,
         water = class11
  ) %>%
  # scale landcover factors here!
  mutate(across(alfalfa:water, ~ as.vector(scale(.))))

landcover8Fixed <- landcoverFixed %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedy = class5,
         wet =  class8,
         water = class11
  ) %>%
  # scale landcover factors here!
  mutate(across(alfalfa:water, ~ as.vector(scale(.))))

# make list of landcover tables to prepare for joins
# define list of names
lcNames <- c('landcover7', 'landcover8', 'landcover7Fixed', 'landcover8Fixed')
lcIn <- list(landcover7, landcover8, landcover7Fixed, landcover8Fixed)
names(lcIn) <- lcNames
# make empty list to hold outputs
landCoverTabs <- list()

# for each version of the landcover tables, join the mean_density_field table.
for (i in 1:length(lcIn)) {

  # Reshape landcover data
  landcover_wide <- lcIn[[i]] %>%
    # finish tidying cols after collapsing klasses
    select(-(class0:class11)) %>%
    # relevel distanceWeights
    mutate(distanceWeight = fct_relevel(distanceWeight,
                                        c('const', 'no'),
                                        after = Inf)) %>%
    # widen
    pivot_wider(id_cols = c(site, field),
                names_from = distanceWeight,
                values_from = alfalfa:water) %>%
    mutate(id = paste0(site, field), .before = site)

  # join landcover data to arthropod density data
  field_data <- left_join(mean_density_field, landcover_wide) %>%
    # tidy: remove extra 'site' and 'field' cols
    select(-site, -field)

  # add to output list
  landCoverTabs[[i]] <- field_data

}

# name list elements
names(landCoverTabs) <- lcNames

# clean up env
rm(data, landcover, landcoverFixed, mean_density_field,
   field_data, landcover7, landcover8, landcover7Fixed, landcover8Fixed)

# define color palette ####
# Acyrthosiphon aphids: #548235
# Non-Acyrthosiphon aphids: #3B3838
# Spring: #76db91
# Fall: #9e3c21
# Predators: 'Spectral'
# Sites: 'Set1'

# Arthropod data summary ####
## Aphid histograms ####
# density plot
PLOT_aphid_density <- mean_density_long %>%
  # focus on aphids only
  filter(Taxa %in% c('Acyrthosiphon', 'NonAcy')) %>%
  # relevel factors
  mutate(Taxa = fct_relevel(Taxa, 'Acyrthosiphon'),
         Season = fct_relevel(Season, 'Spring')) %>%
  ggplot(aes(y = Taxa, x = Mean_Density, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Aphids, Log+1 Transformation',
       x = 'log(Density + 1)',
       y = 'Taxon') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 10) +
  scale_y_discrete(labels=c("Pea aphids +\n Blue alfalfa aphids",
                            "Other taxa")) +
  scale_fill_manual(values = c('#548235', '#3b3838'),
                    guide = 'none')
# plot
PLOT_aphid_density

# boxplot
PLOT_aphid_boxplot <- mean_density_long %>%
  # focus on aphids only
  filter(Taxa %in% c('Acyrthosiphon', 'NonAcy')) %>%
  # relevel factors
  mutate(Taxa = fct_relevel(Taxa, 'Acyrthosiphon'),
         Season = fct_relevel(Season, 'Spring')) %>%
  ggplot(aes(x = Taxa, y = Mean_Density, fill = Taxa)) +
  geom_boxplot() +
  labs(title = 'Aphids, Log+1 Transformation',
       y = 'log(Density + 1)',
       x = 'Taxon') +
  ylim(c(0, 9)) +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, ncol = 2) +
  scale_x_discrete(labels=c("Pea aphids +\n Blue alfalfa aphids",
                            "Other taxa")) +
  scale_fill_manual(values = c('#548235', '#3B3838'),
                    guide = 'none')
# plot
PLOT_aphid_boxplot

# boxplot, all aphids combined, spring vs fall
PLOT_aphidsBySeason_boxplot <- mean_density_long %>%
  # focus on aphids only
  filter(Taxa == 'AllAph') %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, 'Spring')) %>%
  ggplot(aes(x = Season, y = Mean_Density, fill = Season)) +
  geom_boxplot() +
  labs(title = 'Aphids, Log+1 Transformation',
       subtitle = 'Pooled across taxa and split by season',
       x = 'Season',
       y = 'log(Density + 1)') +
  ylim(c(0, 9)) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('#76db91', '#9e3c21'),
                    guide = 'none')
# plot
PLOT_aphidsBySeason_boxplot

## Predator histogram ####
# define list of predators
predlist <- c('Arachnida', 'Anthocoridae', 'Nabis', 'Coccinellidae',
              'Geocoris', 'Ichneumonoidea')
# build density plot
# density plot + boxplot
PLOT_predator_density <- mean_density_long %>%
  # focus on aphids only
  filter(Taxa %in% predlist) %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, 'Spring')) %>%
  ggplot(aes(y = Taxa, x = Mean_Density, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Predators, Log+1 Transformation, all seasons',
       subtitle = 'Plot-level density (means of subplots within a plot)',
       x = 'log(Density + 1)',
       y = 'Taxon') +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Spectral')
# plot
PLOT_predator_density


# Landcover data summary ####
## Alfalfa reclassification effect ####

# make a single large data frame for comparisons
# empty list to hold modified data frames
lcCompares <- list()

for (i in 1:length(lcIn)) {

  lcCompares[[i]] <- lcIn[[i]] %>%
    # lengthen
    pivot_longer(alfalfa:water,
                 names_to = 'klass',
                 values_to = 'weightedArea') %>%
    # drop unused cols
    select(!starts_with('class')) %>%
    # make id col
    unite(id, siteId, distanceWeight, klass, remove = FALSE)

}
# add names
names(lcCompares) <- c('landcover7', 'landcover8',
                 'landcover7Fixed', 'landcover8Fixed')
# join tables with the same number of klasses
sevenClass <- lcCompares$landcover7Fixed %>%
  rename(weightedAreaFixedKlass = weightedArea) %>%
  left_join(lcCompares$landcover7) %>%
  mutate(areaDifference = (weightedArea-weightedAreaFixedKlass)/weightedArea)

eightClass <- lcCompares$landcover8Fixed %>%
  rename(weightedAreaFixedKlass = weightedArea) %>%
  left_join(lcCompares$landcover8)%>%
  mutate(areaDifference = (weightedArea-weightedAreaFixedKlass)/weightedArea)
# put the two versions in a list and add names
klassCompares <- list(sevenClass, eightClass)
names(klassCompares) <- c('sevenClass', 'eightClass')

# plot the change in areas by field
PLOTS_klassCompares <- list()
for (i in 1:length(klassCompares)) {

  # increase in alfalfa due to manual classification
  p <- klassCompares[[i]] %>%
    filter(klass == 'alfalfa') %>%
    mutate(distanceWeight =
             fct_relevel(distanceWeight,
                         c('sig1', 'sig2', 'sig3',
                           'sig4', 'sig5', 'const', 'no')),
           FieldID = reorder_within(siteId,
                                    desc(areaDifference),
                                    distanceWeight,
                                    sep = '_')) %>%
  ggplot(aes(x = FieldID,
               y = areaDifference*-1,
               fill = site)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Increase in alfalfa area score due to manual classification',
         x = 'Field',
         y = 'Area Score change, proportion increase') +
    facet_wrap(~ distanceWeight, scales = 'free') +
    scale_fill_brewer(palette = 'Set1') +
    # scale_y_reordered() +
    scale_x_discrete(label = function(x) str_extract(x, '[:alnum:]+'))

  # decrease in weedy or weedyWet klass due to manual classification
  q <- klassCompares[[i]] %>%
    filter(klass %in% c('weedy', 'weedyWet')) %>%
    mutate(distanceWeight =
             fct_relevel(distanceWeight,
                         c('sig1', 'sig2', 'sig3',
                           'sig4', 'sig5', 'const', 'no')),
           FieldID = reorder_within(siteId,
                                    areaDifference,
                                    distanceWeight,
                                    sep = '_')) %>%
    ggplot(aes(x = FieldID,
               y = areaDifference*-1,
               fill = site)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Decrease in weedy or weedyWet area score
         due to manual classification',
         x = 'Field',
         y = 'Area Score change, proportion decrease') +
    facet_wrap(~ distanceWeight * klass, scales = 'free') +
    scale_fill_brewer(palette = 'Set1') +
    # scale_y_reordered() +
    scale_x_discrete(label = function(x) str_extract(x, '[:alnum:]+'))

  PLOTS_klassCompares[[i]] <- list(p, q)

}

names(PLOTS_klassCompares) <- c('sevenClass', 'eightClass')
# note: alfalfa plots for 7- and 8- class are the same, because the effect of
# moving area TO alfalfa is the same, regardless of what class the area is
# moving FROM
PLOTS_klassCompares$sevenClass[1]
PLOTS_klassCompares$sevenClass[2]
PLOTS_klassCompares$eightClass[2]

# in summary, manual classification seems to change Yerington the most

## weedy/wet binning effect ####
# what is the correlation between weedy and wet cover?
lcIn$landcover8 %>%
  ggplot(aes(x = weedy, y = wet)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_smooth(method = 'lm') +
  facet_wrap(~ distanceWeight, scales = 'free') +
  labs(title = 'Correlation between "weedy" and "wet" classes',
       subtitle = 'across all distance weights') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

# generally poor correlation between these two classes
# bin these because they are biologically similar, not because they are
# correlated

# wet and water correlation
lcIn$landcover8 %>%
  ggplot(aes(x = wet, y = water)) +
  geom_point() +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(~ distanceWeight, scales = 'free') +
  labs(title = 'Correlation between "wet" and "water" classes',
       subtitle = 'across all distance weights') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
# moderate correlation here. this breaks down when using weedyWet (7-class)
lcIn$landcover7 %>%
  ggplot(aes(x = weedyWet, y = water)) +
  geom_point() +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(~ distanceWeight, scales = 'free') +
  labs(title = 'Correlation between "weedyWet" and "water" classes',
       subtitle = 'across all distance weights') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

## pca ####
# define list of distance weights
weightList <- c('sig1', 'sig2', 'sig3', 'sig4', 'sig5', 'const', 'no')
# list to hold PCA outputs
pcaOut <- list()
# for each distance weight, make 4 pcas (one for each dataset)
for (i in 1:length(weightList)) {

  # list to hold filtered dataset
  df <- list()
  # filter each dataset
  for (j in 1:length(lcIn)){

    df[[j]] <- lcIn[[j]] %>%
      filter(distanceWeight == weightList[[i]])

  }
  # list to hold principal components
  comps <- list()
  # 7-class and 8-class need different predictors, so can't for-loop here.
  comps[[1]] <- prcomp(~alfalfa+naturalArid+dirt+ag+impermeable+
                    weedyWet+water, data = df[[1]], scale. = TRUE)
  comps[[2]] <- prcomp(~alfalfa+naturalArid+dirt+ag+impermeable+
                         weedy+wet+water, data = df[[2]], scale. = TRUE)
  comps[[3]] <- prcomp(~alfalfa+naturalArid+dirt+ag+impermeable+
                         weedyWet+water, data = df[[3]], scale. = TRUE)
  comps[[4]] <- prcomp(~alfalfa+naturalArid+dirt+ag+impermeable+
                         weedy+wet+water, data = df[[4]], scale. = TRUE)

  # holds pcas for each dataset
  pcas <- list()
  # list of dataset names for plot titles
  dfNames <- c('7class', '8class',
               '7class, alfalfa fixed', '8class, alfalfa fixed')
  # plot the pca for each dataset
  for(k in 1:length(comps)){

  pcas[[k]] <- autoplot(comps[[k]], data = df[[k]], colour = 'site',
           loadings = TRUE, loadings.label = TRUE) +
    scale_colour_brewer(palette = 'Set1') +
    labs(title = paste0(weightList[[i]], ' ', dfNames[[k]]))

  }
  # save the series of pca plots to the pcaOut list
  pcaOut[[i]] <- pcas

}

# print the the pca plots, four at a time (one distance weight at a time)
for (i in 1:length(pcaOut)){

  do.call('grid.arrange', c(unlist(pcaOut[i], recursive = FALSE), ncol = 2))
  ggsave(paste0('pcaplot', i, '.png'), width = 12, height = 8, units = 'in')

}
# # print all the pca plots at once (hard to read)
# do.call('grid.arrange', c(unlist(pcaOut, recursive = FALSE), ncol = 4))

# Margin data summary ####
## Across sites ####
# Shannon diversity
ggplot(data = vegPlots %>% filter(type == 'Margin'),
       aes(x = site, y = shan, fill = season)) +
  geom_boxplot() +
  labs(title = 'Shannon diversity index',
       x = 'Site', y = 'Diversity') +
  scale_fill_manual(values = c('#76db91', '#9e3c21'))
# Plant species richness
ggplot(data = vegPlots %>% filter(type == 'Margin'),
       aes(x = site, y = rich, fill = season)) +
  geom_boxplot() +
  labs(title = 'Plant species richness',
       x = 'Site', y = 'Richness') +
  scale_fill_manual(values = c('#76db91', '#9e3c21'))
# Total plant cover
ggplot(data = vegPlots %>% filter(type == 'Margin') %>%
         mutate(total_cover = select(., 12:132) %>% rowSums(na.rm = TRUE)),
       aes(x = site, y = total_cover, fill = season)) +
  geom_boxplot() +
  labs(title = 'Plant cover %',
       x = 'Site', y = '%') +
  scale_fill_manual(values = c('#76db91', '#9e3c21'))
### Notes ####
## Not seeing much variation here. May want to shaw by field.
## Remember, Yerington will always be missing. Need to fix colors, factor order,
## etc.

# Model selection ####
## Aphids ~ predators ####
### spring data only ####
#### summary tables ####
# define list of aphid taxa
short_aphlist <- c('Acyrthosiphon',
             'Aphis',
             'Therioaphis', 'NonAcy', 'AllAph')

# create empty list to hold dredges
dredges <- list()
# for each aphid taxon, make a global model and dredge through all predators
for (i in 1:length(short_aphlist)) {

  taxon <- short_aphlist[[i]]
  # note: mean_density_wide already has log-transformed arthropod data
  formula <-  paste0(taxon, '~ Arachnida + Coccinellidae + ',
                     'Ichneumonoidea + Nabis + Geocoris + ',
                     'Anthocoridae + (1|Site:Field)')
  mGlobal <- mean_density_wide %>%
    filter(Season=='Spring') %$%
    lmer(formula(formula),
         na.action = 'na.fail',
         REML = FALSE)

  dredges[[i]] <- dredge(mGlobal) # dredge and store output in list

}
names(dredges) <- short_aphlist

# check data sources
importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(ExpVars = names,
           VarWeights = sw,
           .keep = 'none') %>%
    arrange(ExpVars)
})

names(importance_tab) <- short_aphlist
# make tables of model selection results
springTabs <- lapply(dredges, slice_head, n = 5)
names(springTabs) <- short_aphlist
# show tables
# View(springTabs)

# for each aphid taxon, make variable importance heatmap for all aphids
aphPredVI <- list()
for (i in 1:length(short_aphlist)){

  p <- bind_rows(importance_tab, .id = 'Taxon') %>%
    mutate(VarWeights = as.numeric(VarWeights),
           ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
           Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
    filter(Taxon == short_aphlist[[i]]) %>%
    ggplot(aes(x = reorder(ExpVars, desc(VarWeights)),
               y = VarWeights,
               fill = VarWeights)) +
    geom_col() +
    scale_fill_gradient('',
                        low="blue",
                        high="red",
                        breaks = c(0.3,0.9),
                        labels = c('low','high')) +
    labs(title = paste0(short_aphlist[[i]], ' spring'),
         x = '',
         y = 'Variable importance') +
    theme(legend.position = c(0.95,0.85))

  aphPredVI[[i]] <- p

}
# print all plots
for (i in 1:length(aphPredVI)) {

  print(aphPredVI[[i]])

}
# make a larger heatmap that includes all aphid taxa in one figure

# combine all importance tables into one
aphVI <- rbindlist(importance_tab, idcol = TRUE) %>% as_tibble()
# constuct heatmap
aphPredVISpr <- ggplot(data = aphVI,
                     aes(x = .id, y = ExpVars, fill = VarWeights)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Variable importance, spring',
       subtitle = 'aphid ~ predators + (1|Site:Field), plot-level data') +
  coord_flip()

# average beta coef plots
coefs_tab <- list()
for(i in 1:length(short_aphlist)){
  aphPredCoefs <- dredges[[i]] %>%
    tibble %>%
    select(predlist, weight) %>%
    mutate(across(predlist, ~ multiply_by(.x, weight))) %>%
    summarize(across(.fns = ~ mean(.x, na.rm = TRUE)))
  coefs_tab[[i]] <- aphPredCoefs
}
names(coefs_tab) <- short_aphlist
aphPredCoefs <- rbindlist(coefs_tab, idcol = TRUE) %>%
  as_tibble()

aphPredCoefSpr <- aphPredCoefs %>%
  select(-weight) %>%
  pivot_longer(predlist, names_to = 'term', values_to = 'coeff') %>%
  ggplot(aes(x = .id, y = term, fill = coeff)) +
  geom_tile() +
  geom_text(aes(label = round(coeff, 3))) +
  scale_fill_steps2(low = 'blue', mid = 'white', high = 'red') +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Model-averaged coefficients, spring',
       subtitle = 'aphid ~ predators + (1|Site:Field), plot-level data') +
  coord_flip()

# Correlation plots between aphids and predators
aphPredCorrSpr <- mean_density_wide %>%
  # spring data only
  filter(Season == 'Spring') %>%
  # get only the cols we need
  select(short_aphlist, predlist) %>%
  # create a correlation matrix
  correlate() %>%
  # drop rows and cols to keep only aphid ~ predator correlations
  select(term, short_aphlist) %>%
  filter(term %in% predlist) %>%
  # pivot longer
  pivot_longer(short_aphlist,
               names_to = 'Aphid',
               values_to = 'correlation') %>%
  # build ggplot
  ggplot(aes(x = Aphid, y = term, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2))) +
  scale_fill_steps2(low = 'blue', mid = 'white', high = 'red') +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Correlations, aphids ~ predators, spring',
       subtitle = 'From plot-level data') +
  coord_flip()

# show plots
grid.arrange(aphPredVISpr, aphPredCoefSpr, aphPredCorrSpr, nrow = 3)

##### build models on field-level data for comparison ####
# create empty list to hold dredges
dredges2 <- list()
# for each aphid taxon, make a global model and dredge through all predators
for (i in 1:length(short_aphlist)) {

  taxon <- short_aphlist[[i]]
  # note: mean_density_wide already has log-transformed arthropod data
  formula <-  paste0(taxon, '~ Arachnida + Coccinellidae + ',
                     'Ichneumonoidea + Nabis + Geocoris + ',
                     'Anthocoridae + (1|Site)')
  mGlobal <- landCoverTabs[[1]] %>%
    filter(Season=='Spring') %$%
    lmer(formula(formula),
         na.action = 'na.fail',
         REML = FALSE)

  dredges2[[i]] <- dredge(mGlobal) # dredge and store output in list

}
names(dredges2) <- short_aphlist

# check data sources
importance_tab2 <- lapply(dredges2, function (x) {
  sw(x) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(ExpVars = names,
           VarWeights = sw,
           .keep = 'none') %>%
    arrange(ExpVars)
})

names(importance_tab2) <- short_aphlist
# make tables of model selection results
springTabs2 <- lapply(dredges2, slice_head, n = 5)
names(springTabs2) <- short_aphlist


# combine all importance tables into one
aphVI <- rbindlist(importance_tab2, idcol = TRUE) %>% as_tibble()
# constuct heatmap
aphPredVISpr <- ggplot(data = aphVI,
                       aes(x = .id, y = ExpVars, fill = VarWeights)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Variable importance, spring',
       subtitle = 'aphid ~ predators + (1|Site), field-level data') +
  coord_flip()

# average beta coef plots
coefs_tab <- list()
for(i in 1:length(short_aphlist)){
  aphPredCoefs <- dredges2[[i]] %>%
    tibble %>%
    select(predlist, weight) %>%
    mutate(across(predlist, ~ multiply_by(.x, weight))) %>%
    summarize(across(.fns = ~ mean(.x, na.rm = TRUE)))
  coefs_tab[[i]] <- aphPredCoefs
}
names(coefs_tab) <- short_aphlist
aphPredCoefs <- rbindlist(coefs_tab, idcol = TRUE) %>%
  as_tibble()

aphPredCoefSpr <- aphPredCoefs %>%
  select(-weight) %>%
  pivot_longer(predlist, names_to = 'term', values_to = 'coeff') %>%
  ggplot(aes(x = .id, y = term, fill = coeff)) +
  geom_tile() +
  geom_text(aes(label = round(coeff, 3))) +
  scale_fill_steps2(low = 'blue', mid = 'white', high = 'red') +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Model-averaged coefficients, spring',
       subtitle = 'aphid ~ predators + (1|Site), field-level data') +
  coord_flip()

# Correlation plots between aphids and predators
aphPredCorrSpr <- landCoverTabs[[1]] %>%
  # spring data only
  filter(Season == 'Spring') %>%
  # get only the cols we need
  select(short_aphlist, predlist) %>%
  # create a correlation matrix
  correlate() %>%
  # drop rows and cols to keep only aphid ~ predator correlations
  select(term, short_aphlist) %>%
  filter(term %in% predlist) %>%
  # pivot longer
  pivot_longer(short_aphlist,
               names_to = 'Aphid',
               values_to = 'correlation') %>%
  # build ggplot
  ggplot(aes(x = Aphid, y = term, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2))) +
  scale_fill_steps2(low = 'blue', mid = 'white', high = 'red') +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Correlations, aphids ~ predators, spring',
       subtitle = 'From field-level data') +
  coord_flip()

# show plots
grid.arrange(aphPredVISpr, aphPredCoefSpr, aphPredCorrSpr, nrow = 3)




#
# #### best models ####
# # extract and examine the top Acyrthosiphon model.
# # old now ##
# best.acy.mod <- get.models(springTabs[[1]], subset = 3)[[1]]
# summary(best.acy.mod)
# # must remake to plot effects.
# best.acy.mod.sp <- lmer(Acyrthosiphon ~ Coccinellidae + (1|Site:Field),
#                         data = mean_density_wide %>% filter(Season == 'Spring'),
#                         REML = FALSE, na.action = 'na.omit')
# summary(best.acy.mod.sp)
# # Note: no data seems to be missing here.
#
# # make plot
# # png('spring_acy_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# q <- plot(allEffects(best.acy.mod.sp, residuals = TRUE),
#      main = 'Acyrthosiphon, Spring, nested random effects',
#      id = list(n = 36, labels = mean_density_wide %>%
#                  filter(Season == 'Spring') %>%
#                  unite('id', Site, Field, Plot, sep = '.') %>%
#                  pull(id)))
#
# # dev.off()
#
# ## better plot
# # png('lattice.png', height = 6, width = 6, units = 'in', res = 300)
# trellis.par.set(list(par.xlab.text = list(cex=2),
#                      par.ylab.text = list(cex=2),
#                      par.main.text = list(col = "blue", cex=0.5)))
# plot(effect('Coccinellidae',best.acy.mod.sp,
#             residuals = T),
#      partial.residuals = list(smooth=F),
#      axes = list(x = list(Coccinellidae = list(lab = 'Log(Ladybug density)')),
#                  y = list(lab = 'Log(Aphid density)')),
#      main = NULL,
#      lattice=list(key.args=list(axis.text = list(cex=10))))
# # dev.off()
#
# # # not run:

# # b <- ggemmeans(bmod, terms = 'Coccinellidae')
# #
# # ggplot(b, aes(x, predicted)) +
# #   geom_point()+
# #   geom_line() +
# #   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
#
# # extract and examine the top Aphis model.
# best.aphis.mod <- get.models(springTabs[[2]], subset = 1)[[1]]
# summary(best.aphis.mod)
# # Must remake to plot effects.
# best.aphis.mod <- lmer(Aphis ~ Ichneumonoidea + (1|Site) + (1|Field),
#                        data = mean_density_wide %>% filter(Season == 'Spring'),
#                        REML = FALSE)
# # make plot
# # png('spring_aphis_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# plot(allEffects(best.aphis.mod, residuals = TRUE),
#      main = 'Aphis, Spring',
#      id = list(n = 36, labels = mean_density_wide %>%
#                  filter(Season == 'Spring') %>%
#                  unite('id', Site, Field, Plot, sep = '.') %>%
#                  pull(id)))
# # dev.off()
#
# # extract and examine the top Therioaphis model.
# best.therio.mod <- get.models(springTabs[[3]], subset = 1)[[1]]
# summary(best.therio.mod)
# # Must remake to plot effects.
# best.therio.mod <- lmer(Therioaphis ~ Geocoris + (1|Site) + (1|Field),
#                         data = mean_density_wide %>% filter(Season == 'Spring'),
#                         REML = FALSE)
# # make plot
# # png('spring_therio_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# plot(allEffects(best.therio.mod, residuals = TRUE),
#      main = 'Therioaphis, Spring',
#      id = list(n = 36, labels = mean_density_wide %>%
#                  filter(Season == 'Spring') %>%
#                  unite('id', Site, Field, Plot, sep = '.') %>%
#                  pull(id)))
# # dev.off()

#### plot - level data with all factors ####
# make df of log(total abundance) and add margin and lc data, but don't pool
allFactors_long <- data_long %>%
  # drop 'vial'
  select(-Vial) %>%
  # divide 'Pre-' values by 3 to account for unequal sampling area
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  # log transform
  mutate(Density = log(Density + 1)) %>%
  # pivot wider to make col for each taxon
  pivot_wider(names_from = Taxa, values_from = Density) %>%
  # create unique id for each site:field combo
  mutate(id = paste0(Site, '0', Field)) %>%
  # join margin data
  left_join(field_margins, by = c('id'='id', 'Season'='Season')) %>%
  # add dummy vars for treatment
  cbind(., to.dummy(.$Treatment, 'Trt')) %>%
  # create watering method column
  mutate(wateringMethod = case_when(Site %in% c('Fallon', 'Lovelock') ~
                                      'Flooding',
                                    Site == 'Minden' ~
                                      'Sprinklers',
                                    Site == 'Yerington' & Field %in% c(2, 3) ~
                                      'Flooding',
                                    Site == 'Yerington' & Field == 1 ~
                                      'Sprinklers')) %>%
  # arrange sensibly for readabilicy
  arrange(Season, id, Plot, Treatment) %>%
  # relocate for readability
  relocate(id) %>%
  relocate(starts_with('Trt'), .after = Treatment) %>%
  relocate(wateringMethod, shan, rich, totalCover, .after = Season)

# now, join landcover
allFacts <- list()
for (i in 1:length(landCoverTabs)) {

  lcOneSeason <- landCoverTabs[[i]] %>%
    filter(Season == 'Spring')
  allFacts[[i]] <- left_join(allFactors_long, lcOneSeason, by = 'id') %>%
    # drop duplicate cols
    select(-ends_with('.y')) %>%
    # strip '.x' suffix
    rename_with(~ str_sub(.x, end = -3), .cols = ends_with('.x')) %>%
    tibble

}
names(allFacts) <- names(landCoverTabs)

#### review sham trt effects ####
# overall effect on predators
# make sham-cont column
# this is inappropriate use of logs


# try this again with unlogged data
diffData_wide <- data_long %>%
  filter(Treatment != 'Pre-') %>%
  select(-Vial) %>%
  pivot_wider(names_from = Taxa, values_from = Density) %>%
  pivot_wider(names_from = Treatment, values_from = Arachnida:NonAcy) %>%
  mutate(Arachnida = Arachnida_Sham - Arachnida_Control,
         Nabis = Nabis_Sham - Nabis_Control,
         Geocoris = Geocoris_Sham - Geocoris_Control,
         Anthocoridae = Anthocoridae_Sham - Anthocoridae_Control,
         Coccinellidae = Coccinellidae_Sham - Coccinellidae_Control,
         Ichneumonoidea = Ichneumonoidea_Sham - Ichneumonoidea_Control,
         Acyrthosiphon = Acyrthosiphon_Sham - Acyrthosiphon_Control,
         Aphis = Aphis_Sham - Aphis_Control,
         Therioaphis = Therioaphis_Sham - Therioaphis_Control,
         Lygus = Therioaphis_Sham - Therioaphis_Control,
         Thysanoptera = Thysanoptera_Sham - Thysanoptera_Control,
         Other = Other_Sham - Other_Control,
         AllAph = AllAph_Sham - AllAph_Control,
         NonAcy = NonAcy_Sham - NonAcy_Control) %>%
  select(-ends_with('_Sham'), -ends_with('_Control'))

diffData <- diffData_wide %>%
  pivot_longer(cols = Arachnida:NonAcy,
               names_to = 'Taxon',
               values_to = 'Difference')

##### predators #####
predSpringStats <- list()
for(i in 1:length(predlist)){

  pred <- as.name(predlist[[i]])
  t <- t.test(diffData_wide %>%
                filter(Season == 'Spring') %>%
                select(all_of(pred)))
  stats <- c(t$estimate, t$conf.int[[1]], t$conf.int[[2]], t$p.value)
  predSpringStats[[i]] <- stats

}
names(predSpringStats) <- predlist

predSpringStatsDf <- as.data.frame(do.call(rbind, predSpringStats)) %>%
  rename(mean = `mean of x`, lower = V2, upper = V3, p = V4) %>%
  rownames_to_column('taxon')

predFallStats <- list()
for(i in 1:length(predlist)){

  pred <- as.name(predlist[[i]])
  t <- t.test(diffData_wide %>%
                filter(Season =='Fall') %>%
                select(all_of(pred)))
  stats <- c(t$estimate, t$conf.int[[1]], t$conf.int[[2]], t$p.value)
  predFallStats[[i]] <- stats

}
names(predFallStats) <- predlist

predFallStatsDf <- as.data.frame(do.call(rbind, predFallStats)) %>%
  rename(mean = `mean of x`, lower = V2, upper = V3, p = V4) %>%
  rownames_to_column('taxon')

predStats2 <- list()
for (i in 1:length(predlist)) {
  t<-lm(eval(as.name(predlist[[i]])) ~ Season, data = diffData_wide)
  summary <- summary(t)
  stats <- c(t$coefficients[[2]], summary$coefficients[[2,4]])
  predStats2[[i]] <- stats
  }
names(predStats2) <- predlist
predStatsDf2 <- as.data.frame(do.call(rbind, predStats2)) %>%
  rename(FallEffect = V1, p = V2) %>%
  rownames_to_column('taxon')

# plot. uses predSpringStatsDf (or fall version) for mean and error
diffData %>%
  filter(Taxon %in% predlist,
         Season == 'Spring') %>%
  ggplot(aes(x = Taxon, y = Difference, fill = Taxon)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Spectral') +
  geom_point(inherit.aes = FALSE,
             data = predSpringStatsDf,
             aes(x = taxon, y = mean),
             position = position_nudge(x = 0.2)) +
  geom_errorbar(inherit.aes = FALSE,
                data = predSpringStatsDf,
                position = position_nudge(x = 0.2),
    aes(ymin = lower, ymax = upper, x = taxon),
    width = 0.2) +
  geom_hline(yintercept = 0, color = 'red') +
  theme(legend.position = 'none') +
  coord_flip(ylim = c(-10, 14))
  #

predStatsDf2 %>%
  # select(taxon, mean, lower, upper, p) %>%
  arrange(desc(taxon)) %>%
  mutate(across(where(is.double), ~round(.x, 3))) %>%
  grid.table()

predFallStatsDf %>%
  select(taxon, mean, lower, upper, p) %>%
  arrange(desc(taxon)) %>%
  mutate(across(where(is.double), ~round(.x, 3))) %>%
  grid.table()

diffData %>%
  filter(Season == 'Fall',
         Taxon == 'Ichneumonoidea') %>%
  ggplot(aes(x = Taxon, y = Difference)) +
  geom_boxplot(fill ='#99d594') +
  geom_point(inherit.aes = FALSE,
             data = predFallStatsDf %>% filter(taxon == 'Ichneumonoidea'),
             aes(x = taxon, y = mean),
             position = position_nudge(x = 0.2)) +
  geom_errorbar(inherit.aes = FALSE,
                data = predFallStatsDf %>% filter(taxon == 'Ichneumonoidea'),
                position = position_nudge(x = 0.2),
                aes(ymin = lower, ymax = upper, x = taxon),
                width = 0.2) +
  geom_hline(yintercept = 0, color = 'red') +
  theme(legend.position = 'none') +
  coord_flip()


##### aphids #####

short_aphlist <- c('Acyrthosiphon', 'Aphis', 'Therioaphis', 'NonAcy', 'AllAph')

aphSpringStats <- list()
for(i in 1:length(short_aphlist)){

  aph <- as.name(short_aphlist[[i]])
  t <- t.test(diffData_wide %>%
                filter(Season =='Spring') %>%
                select(all_of(aph)))
  stats <- c(t$estimate, t$conf.int[[1]], t$conf.int[[2]], t$p.value)
  aphSpringStats[[i]] <- stats

}
names(aphSpringStats) <- short_aphlist
View(diffData_wide)
aphSpringStatsDf <- as.data.frame(do.call(rbind, aphSpringStats)) %>%
  rename(mean = `mean of x`, lower = V2, upper = V3, p = V4) %>%
  rownames_to_column('taxon')

aphFallStats <- list()
for(i in 1:length(short_aphlist)){

  pred <- as.name(short_aphlist[[i]])
  t <- t.test(diffData_wide %>%
                filter(Season =='Fall') %>%
                select(all_of(pred)))
  stats <- c(t$estimate, t$conf.int[[1]], t$conf.int[[2]], t$p.value)
  aphFallStats[[i]] <- stats

}
names(aphFallStats) <- short_aphlist

aphFallStatsDf <- as.data.frame(do.call(rbind, aphFallStats)) %>%
  rename(mean = `mean of x`, lower = V2, upper = V3, p = V4) %>%
  rownames_to_column('taxon')



# plot. uses predSpringStatsDf (or fall version) for mean and error
diffData %>%
  filter(Season == 'Fall',
         Taxon %in% short_aphlist) %>%
  ggplot(aes(x = Taxon, y = Difference, fill = Taxon)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Spectral') +
  geom_point(inherit.aes = FALSE,
             data = aphFallStatsDf,
             aes(x = taxon, y = mean),
             position = position_nudge(x = 0.2)) +
  geom_errorbar(inherit.aes = FALSE,
                data = aphFallStatsDf,
                position = position_nudge(x = 0.2),
                aes(ymin = lower, ymax = upper, x = taxon),
                width = 0.2) +
  geom_hline(yintercept = 0, color = 'red') +
  theme(legend.position = 'none') +
  coord_flip()

#
View(diffData)
aphSpringStatsDf %>%
  select(taxon, mean, lower, upper, p) %>%
  arrange(desc(taxon)) %>%
  mutate(across(where(is.double), ~round(.x, 3))) %>%
  grid.table()

aphFallStatsDf %>%
  select(taxon, mean, lower, upper, p) %>%
  arrange(desc(taxon)) %>%
  mutate(across(where(is.double), ~round(.x, 3))) %>%
  grid.table()

##### correlation of aphids and predators #####
diffData_wide %>%
  filter(AllAph > -3000) %>%
  pivot_longer(predlist, names_to = 'Predator', values_to = 'Difference') %>%
  ggplot(aes(x = AllAph, y = Difference, fill = Predator, color = Predator)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_fill_brewer(palette = 'Spectral') +
  scale_color_brewer(palette = 'Spectral') +
  labs(title = 'All Seasons',
       subtitle = 'Correlations between changes in aphid density in sham treatments
       and changes in predator density in sham treatments')



# Correlation plots between aphids and predators
diffData_wide %>%
  # spring data only
  filter(Season == 'Spring') %>%
  # get only the cols we need
  select(predlist) %>%
  # create a correlation matrix
  correlate() %>%
  shave %>%
  rplot(print_cor = TRUE)

diffData_wide %>%
  filter(Season == 'Spring') %>%
  ggplot(aes(Anthocoridae, Coccinellidae)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_poly_eq()

diffData %>%
  filter(Season == 'Fall',
         Taxon %in% predlist) %>%
  # create watering method column
  mutate(wateringMethod = case_when(Site %in% c('Fallon', 'Lovelock') ~
                                      'Flooding',
                                    Site == 'Minden' ~
                                      'Sprinklers',
                                    Site == 'Yerington' & Field %in% c(2, 3) ~
                                      'Flooding',
                                    Site == 'Yerington' & Field == 1 ~
                                      'Sprinklers')) %>%
ggplot(aes(x = Taxon, y = Difference, fill = Taxon)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Spectral') +
  geom_hline(yintercept = 0, color = 'red') +
  theme(legend.position = 'none') +
  facet_wrap(~wateringMethod) +
  coord_flip(ylim = c(-10, 14))


allFactors_long %>%
  filter(Season == 'Spring') %>%
  pivot_longer(predlist, names_to = 'Predator', values_to = 'Density') %>%
  ggplot(aes(AllAph, Density, fill = Predator, color = Predator)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_poly_eq() +
  labs(title = 'All sites and seasons, subplot-level data, log-transformed')



#### Beths Q ####
# compare intercept in these models, per Beth's advice
# not sure whether mixed mod is appropriate - must resolve this
##### first with log-transformed data ####
qData <- allFactors_long %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-')
# # remove outlier aphid count
# filter(AllAph < 8.2936) ## DON'T DO THIS!!
# MUST REMOVE ENTIRE PLOT
# find offending plot
qData %>% filter(AllAph >= 8.2936)
nrow(qData)
# Lovelock 2-2
qData %<>% filter(!(Site == "Lovelock" & Field == 2 & Plot == 2))
nrow(qData)

# filter where sham was 0
qNonZero <- qData %>%
  filter(!(Coccinellidae == 0 & Treatment == 'Sham'))

qZero <- qData %>%
  filter((Coccinellidae == 0 & Treatment == 'Sham'))
View(diffData_wide %>% filter(Season == "Spring"))


# make lm and lmer mods
linMod <- lm(AllAph ~ Coccinellidae + Treatment, data = qData)
linModNonZero <- lm(AllAph ~ Coccinellidae + Treatment, data = qNonZero)
linModZero <- t.test(qZero$AllAph)

linModMixed <- lmer(AllAph ~ Coccinellidae + Treatment + (1|Site:Field),
                  data = qData)
# plot effs
plot(allEffects(linMod))
plot(allEffects(linModNonZero))
plot(allEffects(linModMixed))

# need ggiraphExtra for the following

# this seems to work well for the non-mixed mod
ggAncova(linMod, interactive = TRUE)
ggAncova(linModNonZero, interactive = TRUE)
ggAncova(linModZero, interactive = TRUE)


tab_model(linMod)
tab_model(linModNonZero)
tab_model(linModMixed)

# do models on jenks categories
# need sumQData from below
qJoin <- qData %>%
  left_join(sumQData)

linModLow <- qJoin %>%
  filter(logDensJenks == 'Low') %$%
  lm(AllAph ~ Coccinellidae + Treatment)

linModHigh <- qJoin %>%
  filter(logDensJenks == 'High') %$%
  lm(AllAph ~ Coccinellidae + Treatment)

ggAncova(linModLow)
ggAncova(linModHigh)

tab_model(linModLow)
tab_model(linModHigh)

# This is beth's fave approach: ####
# do models on +/- sham effect
qDiffJoin <- qJoin %>%
  left_join(diffData_wide %>%
              filter(Season == 'Spring') %>%
              mutate(diffLadybugs = Coccinellidae) %>%
              select(Site, Field, Plot, diffLadybugs)) %>%
  mutate(shamEffect = case_when(diffLadybugs > 0 ~ "positive",
                                diffLadybugs <= 0 ~ "neutral/negative"))


linModPos <- qDiffJoin %>%
  filter(shamEffect == 'positive') %$%
  lm(AllAph ~ Coccinellidae + Treatment)

linModNeg <- qDiffJoin %>%
  filter(shamEffect == 'neutral/negative') %$%
  lm(AllAph ~ Coccinellidae + Treatment)

ggAncova(linModPos)
ggAncova(linModNeg)

tab_model(linModPos)
tab_model(linModNeg)

# looking good
# repeat for Anthocoridae, Arachnida
# repeat across fall, sp+fa combined

diffMod <- lm(AllAph ~ Coccinellidae, data = diffData_wide %>%
                filter(Season == 'Spring'))
tab_model(diffMod)
##### developing composite predator effect #####
# could sum all predators
qData %>% mutate(AllPreds = Anthocoridae + Arachnida + Coccinellidae) %>%
  ggplot(aes(x = AllPreds, y = AllAph)) +
  geom_point(aes(color = Treatment)) +
  geom_smooth(method = 'lm') +
  labs(title = '3 preds summed')

qData %>% mutate(AllPreds = Anthocoridae + Arachnida + Coccinellidae) %$%
  lm(AllAph ~ AllPreds + Treatment) -> allPredMod

ggAncova(allPredMod)

tab_model(allPredMod)

# could scale and then sum
qData %>% mutate(across(.cols = c(Anthocoridae, Arachnida, Coccinellidae),
                        ~ scale(.x, center = FALSE),
                        .names = '{.col}_scale')) %>%
  mutate(AllPreds = Anthocoridae_scale + Arachnida_scale + Coccinellidae_scale) %>%
  ggplot(aes(x = AllPreds, y = AllAph)) +
  geom_point(aes(color = Treatment)) +
  geom_smooth(method = 'lm') +
  labs(title = 'preds scaled then summed')



# is sham effect similar across predators where it is known to have an effect



# fixed int: THIS IS NONSENSE
# You could subtract the explicit intercept from the regressand and then fit the
# intercept-free model:

intercept <- linMod$coefficients[[1]]
linModFix <- lm(I(AllAph - intercept) ~
                  0+Coccinellidae * Treatment, data = qData)

ggAncova(linModFix)

ggplot(qData, aes(y = AllAph, x = Coccinellidae, color = Treatment))+
  geom_point()+
  geom_smooth(method = 'lm', formula = y~0+x, fullrange = T)+
  xlim(0,3)

tab_model(linModFix)
## END NONSENSE

# ggAncova(linModMixed, interactive = TRUE) # not good

Anova(linMod)
summary(linMod)

tab_model(linMod)
tab_model(linMod, file = 'linMod.html')
webshot('linMod.html', 'linMod.png')


tryAov <- aov(AllAph ~ Coccinellidae + Treatment, data = qData)
Anova(tryAov)

shamDat <- qData %>% filter(Treatment == 'Sham')
contDat <- qData %>% filter(Treatment == 'Control')

nrow(contDat)

tryAov2 <- aov(AllAph ~ Treatment, data = qData)
Anova(tryAov2)


# this is the same as the lm, as predicted
# drop coccinellidae?
tryAov <- aov(AllAph ~ Treatment, data = qData)
Anova(tryAov)
summary(tryAov)


wideQ <- qData %>%
  select(-`Trt.Pre-`, -`Trt.Control`, -`Trt.Sham`) %>%
  pivot_wider(names_from = Treatment, values_from = shan:NonAcy)

tryT <- t.test(wideQ$AllAph_Control, wideQ$AllAph_Sham)
tryTPaired <- t.test(wideQ$AllAph_Control, wideQ$AllAph_Sham, paired = TRUE)


ggplot(qData, aes(x =Treatment, y=AllAph)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))

# try interaction - sham can't attract ladybugs if they aren't present
# restore outlier
qData <- allFactors_long %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-')

interLm <- lm(AllAph ~ Coccinellidae * Treatment, data = qData)
summary(interLm)
plot(allEffects(interLm))
ggEffect(AllAph ~ Coccinellidae * Treatment, data = qData)


#### Jenks #####
# given log-transform && keeping the outlier, what about thresholds of ladybug
# density?
# restore outlier into data
qData <- allFactors_long %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-')

# find jenks or kmeans clusters
install.packages('BAMMtools')

install.packages('classInt')


# summarize data to get plot-level means
sumQData <- qData %>%
  group_by(Site, Field, Plot) %>%
  summarize(meanLadybugs = mean(Coccinellidae))
sumRData <- rData %>%
  group_by(Site, Field, Plot, Season) %>%
  summarize(meanLadybugs = mean(Coccinellidae))

# jenks
logJenks <- getJenksBreaks(sumQData$meanLadybugs, 3)
Jenks <- getJenksBreaks(sumRData$meanLadybugs, 3)
# kmeans
logKMeans <- classIntervals(sumQData$meanLadybugs, 2, style = "kmeans")
KMeans <- classIntervals(sumRData$meanLadybugs, 2, style = "kmeans")
# update data tables with category
sumQData %<>%
  mutate(logDensJenks = case_when(meanLadybugs <= logJenks[[1]] ~ 'Zero',
                               meanLadybugs > 0 & meanLadybugs <= logJenks[[2]]
                               ~ 'Low',
                               meanLadybugs > logJenks[[2]] ~ 'High')) %>%
mutate(logDensKMeans = case_when(meanLadybugs <= logKMeans$brks[[1]] ~ 'Zero',
                             meanLadybugs > 0 & meanLadybugs <= logKMeans$brks[[2]]
                             ~ 'Low',
                             meanLadybugs > logKMeans$brks[[2]] ~ 'High'))

sumRData %<>%
  mutate(densJenks = case_when(meanLadybugs <= Jenks[[1]] ~ 'Zero',
                               meanLadybugs > 0 & meanLadybugs <= Jenks[[2]]
                               ~ 'Low',
                               meanLadybugs > Jenks[[2]] ~ 'High'))%>%
  mutate(densKMeans = case_when(meanLadybugs <= KMeans$brks[[1]] ~ 'Zero',
                                meanLadybugs > 0 & meanLadybugs <= KMeans$brks[[2]]
                                ~ 'Low',
                                meanLadybugs > KMeans$brks[[2]] ~ 'High'))

# join to original dfs





# nice plot for logged data
grob <- grobTree(textGrob("K Means breaks", x=0.1,  y=0.9, hjust=0,
                          gp=gpar(col="red", fontsize=13, fontface="italic")))
grob2 <- grobTree(textGrob("Jenks breaks", x=0.9,  y=0.9, hjust=1,
                           gp=gpar(col="blue", fontsize=13, fontface="italic")))
ggplot(sumQData, aes('', meanLadybugs)) +
  geom_violin() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1,
               binwidth = max(sumQData$meanLadybugs)/20,
               aes(fill = logDensJenks, group = logDensJenks)) +
  geom_hline(yintercept = c(logKMeans$brks[[1]] + 0.065, # manual jitter
                            logKMeans$brks[[2]] + 0.01,
                            logKMeans$brks[[3]]+ 0.02), color = 'red') +
  geom_hline(yintercept = c(logJenks[[1]] + 0.075,
                            logJenks[[2]] + 0.075,
                            logJenks[[3]] + 0.03), color = 'blue') +
  annotation_custom(grob) +
  annotation_custom(grob2) +
  labs(y = 'log(Coccinellidae + 1) Density', x = '',
       fill = 'Jenks group',
       title = 'Plot grouping by log(ladybug density)',
       subtitle = 'Spring only') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# nice plot for unlogged data
ggplot(sumRData, aes('', meanLadybugs)) +
  geom_violin() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1,
               binwidth = max(sumRData$meanLadybugs)/20,
               aes(fill = densJenks, group = densJenks)) +
  geom_hline(yintercept = c(KMeans$brks[[1]] + 0.22, # manual jitter
                            KMeans$brks[[2]] - 0.03,
                            KMeans$brks[[3]]+ 0.22), color = 'red') +
  geom_hline(yintercept = c(Jenks[[1]] + 0.27,
                            Jenks[[2]] + 0.27,
                            Jenks[[3]] + 0.27), color = 'blue') +
  annotation_custom(grob) +
  annotation_custom(grob2) +
  labs(y = 'Coccinellidae Density', x = '',
       fill = 'Jenks group',
       title = 'Plot grouping by ladybug density - no log transform') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


# does sham have a larger effect when ladybugs are nonzero?
diffData %>%
  filter(Taxon == 'Coccinellidae',
         Season == 'Spring') %>%
  left_join(sumQData %>% select(Site, Field, Plot, meanLadybugs, logDensJenks)) %>%
  ggplot(aes(x = logDensJenks, y = Difference, fill = logDensJenks)) +
  geom_boxplot() +
  # geom_point(inherit.aes = FALSE,
  #            data = predSpringStatsDf,
  #            aes(x = taxon, y = mean),
  #            position = position_nudge(x = 0.2)) +
  # geom_errorbar(inherit.aes = FALSE,
  #               data = predSpringStatsDf,
  #               position = position_nudge(x = 0.2),
  #               aes(ymin = lower, ymax = upper, x = taxon),
  #               width = 0.2) +
  geom_hline(yintercept = 0, color = 'red') +
  coord_flip(ylim = c(-10, 14))

diffData %>%
  filter(Taxon == 'Coccinellidae',
         Season == 'Spring')  %>%
  left_join(sumQData %>% select(Site, Field, Plot, meanLadybugs, logDensJenks)) ->t#
t%>% filter(Season == 'Spring') %>%
  ggplot(aes(x = meanLadybugs, y = Difference)) +
  geom_jitter(height = 0, width = 0.05, shape = 21) +
  geom_smooth(method = 'lm', formula=y~x+0) +
  labs(x = "Plot-level mean ladybug density (jittered)",
       y = 'Sham - Control ladybug density')

diffMod <- diffData %>%
  filter(Taxon == 'Coccinellidae') %>%
  left_join(sumRData %>% select(Site, Field, Plot, Season, meanLadybugs, densJenks)) %>%
  filter(Season == 'Spring')%$%
  lm(Difference ~ 0 + meanLadybugs)

summary(diffMod)
tab_model(diffMod)

plot(diffMod)


 ##### not log-transformed data ####
allFactors_long_noLog <- data_long %>%
  # drop 'vial'
  select(-Vial) %>%
  # divide 'Pre-' values by 3 to account for unequal sampling area
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  # log transform
  # mutate(Density = log(Density + 1)) %>%
  # pivot wider to make col for each taxon
  pivot_wider(names_from = Taxa, values_from = Density) %>%
  # create unique id for each site:field combo
  mutate(id = paste0(Site, '0', Field)) %>%
  # join margin data
  left_join(field_margins, by = c('id'='id', 'Season'='Season')) %>%
  # add dummy vars for treatment
  cbind(., to.dummy(.$Treatment, 'Trt')) %>%
  # create watering method column
  mutate(wateringMethod = case_when(Site %in% c('Fallon', 'Lovelock') ~
                                      'Flooding',
                                    Site == 'Minden' ~
                                      'Sprinklers',
                                    Site == 'Yerington' & Field %in% c(2, 3) ~
                                      'Flooding',
                                    Site == 'Yerington' & Field == 1 ~
                                      'Sprinklers')) %>%
  # arrange sensibly for readabilicy
  arrange(Season, id, Plot, Treatment) %>%
  # relocate for readability
  relocate(id) %>%
  relocate(starts_with('Trt'), .after = Treatment) %>%
  relocate(wateringMethod, shan, rich, totalCover, .after = Season)

# prep data
rData <- allFactors_long_noLog %>%
  filter(Treatment != 'Pre-')
nrow(rData)

# MUST REMOVE ENTIRE PLOT
# Lovelock 2-2
rData %<>% filter(!(Site == "Lovelock" & Field == 2 & Plot == 2))
pData <- rData %>% filter((Site == "Lovelock" & Field == 2 & Plot == 2))
nrow(rData)
# make mods
linMod_noLog <- lm(AllAph ~ Coccinellidae + Treatment, data = rData)
linModMixed <- lmer(AllAph ~ Coccinellidae + Treatment + (1|Site:Field),
                    data = rData)

ggAncova(linMod_noLog, interactive = TRUE)
ggEffect(linMod_noLog, interactive = TRUE)
ggEffect(linMod, interactive = TRUE)


tab_model(linMod_noLog, file = 'linmodNoLog.html')
webshot('linmodNoLog.html', 'linModNoLog.png')
summary(linMod)

anova(linMod)
Anova(linMod)
Anova(linMod, type = "III")

t <- allFacts[[1]] %>%
  filter(Season == 'Spring')

ichMod <- lmer(AllAph ~ Ichneumonoidea + (1|Site:Field), data = t, REML = FALSE)
cocMod <- lmer(AllAph ~ Coccinellidae + (1|Site:Field), data = t, REML = FALSE)

summary(ichMod)
summary(cocMod)

plot(allEffects(ichMod))
plot(allEffects(cocMod))
allEffects(ichMod)
str(t)

noCovMod <- lm(AllAph ~ Coccinellidae, data = qData)
ggAnova(noCovMod)

'meanDens' <- mean_density_long %>%
  # focus on aphids only
  filter(Taxa %in% predlist) %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, 'Spring')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density)

meanDensTop <- list()
meanDensBottom <- list()
meanDensMiddle <- list()
for (i in 1:length(predlist)) {

  df <- meanDens %>% arrange(eval(as.name(predlist[[i]])))
  meanDensTop[[i]] <-df %>% slice_head(n = 21) %>% select(Site,
                                                          Field,
                                                          Plot,
                                                          Season)
  meanDensBottom[[i]]<- df %>% slice_tail(n = 21) %>% select(Site,
                                                             Field,
                                                             Plot,
                                                             Season)
  meanDensMiddle[[i]]<-df %>% filter(between(row_number(), 22, 44)) %>% select(Site,
                                                            Field,
                                                            Plot,
                                                            Season)

}


names(meanDensTop) <- predlist
names(meanDensBottom) <- predlist

mean_density_wide  %>%
  summarise(Acyrthosiphon = mean(Acyrthosiphon))


lowList <- list()
for (i in 1:length(predlist)){
p<-diffData_wide %>%
  right_join(meanDensTop[[i]], by = c('Site', 'Field', 'Plot', 'Season')) %>%
  filter(AllAph > -3000) %>%
  ggplot(aes(x = eval(as.name(predlist[[i]])))) +
  geom_boxplot() +
  labs(x = predlist[[i]],
       title = paste('Sham effect in the LOWEST tertile of density:',
                     predlist[[i]]),
       subtitle = 'all seasons')# +
  # geom_smooth(method = 'lm') +
  # scale_fill_brewer(palette = 'Spectral') +
  # scale_color_brewer(palette = 'Spectral') +
  # labs(title = paste('All Seasons', predlist[[i]]),
  #      subtitle = 'Correlations between changes in aphid density in sham treatments
  #      and changes in predator density in sham treatments, HIGHEST TERTILE')
print(p)
lowList[[i]]<-p
}

highList <- list()
for (i in 1:length(predlist)){
  p<-diffData_wide %>%
    right_join(meanDensBottom[[i]], by = c('Site', 'Field', 'Plot', 'Season')) %>%
    filter(AllAph > -3000) %>%
    ggplot(aes(x = eval(as.name(predlist[[i]])))) +
    geom_boxplot() +
    labs(x = predlist[[i]],
         title = paste('Sham effect in the HIGHEST tertile of density:',
                       predlist[[i]]),
         subtitle = 'all seasons')# +
  # geom_smooth(method = 'lm') +
  # scale_fill_brewer(palette = 'Spectral') +
  # scale_color_brewer(palette = 'Spectral') +
  # labs(title = paste('All Seasons', predlist[[i]]),
  #      subtitle = 'Correlations between changes in aphid density in sham treatments
  #      and changes in predator density in sham treatments, HIGHEST TERTILE')
  print(p)
highList[[i]]<-p
}

midList <- list()
for (i in 1:length(predlist)){
  p<-diffData_wide %>%
    right_join(meanDensMiddle[[i]], by = c('Site', 'Field', 'Plot', 'Season')) %>%
    filter(AllAph > -3000) %>%
    ggplot(aes(x = eval(as.name(predlist[[i]])))) +
    geom_boxplot() +
    labs(x = predlist[[i]],
         title = paste('Sham effect in the MIDDLE tertile of density:',
                       predlist[[i]]),
         subtitle = 'all seasons')# +
  # geom_smooth(method = 'lm') +
  # scale_fill_brewer(palette = 'Spectral') +
  # scale_color_brewer(palette = 'Spectral') +
  # labs(title = paste('All Seasons', predlist[[i]]),
  #      subtitle = 'Correlations between changes in aphid density in sham treatments
  #      and changes in predator density in sham treatments, HIGHEST TERTILE')
  print(p)
  midList[[i]]<-p
}

for(i in 1:6){

  print(grid.arrange(lowList[[i]], midList[[i]], highList[[i]]))

}

diffData_wide %>%
  right_join(meanDensBottom$Anthocoridae,
             by = c('Site', 'Field', 'Plot', 'Season')) %>%
  arrange(desc(Anthocoridae))


#### try some predator mods that include trt ####
##### now, find the best landcover + trt + aphid model for each pred!! ####
buildLandcoverModTab <- function(taxon = 'empty', data = 'empty',
                                 dataSet = 'empty', m.max = 3){
  # taxon='Arachnida'
  # dataSet = allFacts[1]
  # data= allFacts[[2]]
  # m.max=3

  distList <- c('_no ',
                '_const ',
                '_sig1 ',
                '_sig2 ',
                '_sig3 ',
                '_sig4 ',
                '_sig5 ')

  varList <- data %>%
    select(contains('_')) %>%
    names() %>% str_extract('[:alpha:]+') %>%
    unique()



  if (taxon == 'empty' | !is_tibble(data)) {
    stop(red("Please specify taxon and data \n"), call. = FALSE)
  }

  # print inputs
  cat(yellow('Taxon:'),
      green(taxon),
      yellow('Data:'),
      green(names(dataSet)),
      '\n')
  # get dfname
  dfname <- as.name(deparse(substitute(data)))

  aphList <- append(short_aphlist, 'AllAph')
  cand_mod_tabs_aph <- list()
  for (i in 1:length(aphList)) {
    # i=1
    message(blue('fitting', aphList[[i]], 'models'))
    cand_mod_tabs_dist <- list()
    for (j in 1:length(distList)) {
      # j=1
      # incase full model fails to fit, try 'tricking' dredge
      message(blue('fitting', distList[[j]],'models'))
      # fit reduced model
      fmod.red <- lmer(as.formula(
        paste0(taxon, ' ~ (1|Site:Field)')),
        data = data,
        REML = FALSE,
        na.action = 'na.fail')
      # define full model formula
      form <-formula(
        paste0(taxon, ' ~ ',
               paste0(varList, distList[[j]], '+ ', collapse = ''),
               aphList[[i]],'+ (1|Site:Field)'
        ))
      # Replace reduced model formula with full global model formula.
      attr(fmod.red@frame, "formula") <- form
      # Run dredge() with m.max parameter to avoid convergence failures.
      # fmod.red@call$data <- dfname

      cand_mod_tabs_dist[[j]] <-  # superassign?
        model.sel(lapply(
          dredge(fmod.red, m.lim = c(NA, m.max),
                 trace = 2, evaluate = FALSE),
          eval))

      message(blue(nrow(cand_mod_tabs_dist[[j]]), 'models fit\n'))

    }
    # Rbind the elements of the list together. This forces recalculation of AICc
    cand_mod_tabs_aph[[i]] <- rbind(cand_mod_tabs_dist[[1]],
                                    cand_mod_tabs_dist[[2]],
                                    cand_mod_tabs_dist[[3]],
                                    cand_mod_tabs_dist[[4]],
                                    cand_mod_tabs_dist[[5]],
                                    cand_mod_tabs_dist[[6]],
                                    cand_mod_tabs_dist[[7]])

  }
  full_mod_table <- rbind(cand_mod_tabs_aph[[1]],
                          cand_mod_tabs_aph[[2]],
                          cand_mod_tabs_aph[[3]],
                          cand_mod_tabs_aph[[4]],
                          cand_mod_tabs_aph[[5]])
  return(full_mod_table)
}

# # make full mod tabs for each taxon and version of the data (heavy!)
# predModTabsSpring <- list()
# for (i in 1:length(allFacts)) {
#   message(red('LC Version:', names(allFacts[i])))
#   data <- allFacts[[i]] %>% filter(Season == 'Spring')
#   lcVers <- list()
#   for (j in 1:length(predlist)) {
#     message(red('LC Version:', names(allFacts[i])),
#             blue('Response:', predlist[[j]]))
#     lcVers[[j]] <- buildLandcoverModTab(predlist[[j]], data)
#
#   }
#   names(lcVers) <- names(predlist)
#   predModTabsSpring[[i]] <- lcVers
#
# }
#
# names(predModTabsSpring) <- names(allFacts)
#
# names(predModTabsSpring$landcover7) <- predlist
# names(predModTabsSpring$landcover8) <- predlist
# names(predModTabsSpring$landcover7Fixed) <- predlist
# names(predModTabsSpring$landcover8Fixed) <- predlist
#
# saveRDS(predModTabsSpring, 'predModTabsSpring')

predModTabsSpring <- readRDS('predModTabsSpring')


importance_tab <- sw(predModTabsSpring$landcover8$Coccinellidae) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = "_") %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `sig1` = 'Very aggressive',
                                       `sig2` = 'Aggressive',
                                       `sig3` = 'Moderately aggressive',
                                       `sig4` = 'Moderate',
                                       `sig5` = 'Slight',
                                       `sig6` = 'Minimal',
                                       `const` = 'Constant',
                                       `no` = 'None'))) %>%
  mutate(distWeight = fct_relevel(distWeight, 'Constant', 'None', after = Inf),
         class = recode_factor(class,
                        ag = 'Agricultural',
                        alfalfa = 'Alfalfa',
                        dirt = 'Bare soil +\n dirt road',
                        impermeable = 'Impermeable',
                        naturalArid = 'Natural',
                        water = 'Surface\nwater',
                        weedy = 'Weedy',
                        wet = 'Riparian',
                        AllAph = 'All aphids',
                        Acyrthosiphon = 'Acyrthosiphon',
                        NonAcy = 'Non-Acyrthosiphon',
                        Aphis = 'Aphis',
                        Therioaphis = 'Therioaphis',
                        .ordered = TRUE))

ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(angle = 45))+
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       title = paste0('log(', ')',
                      ' Variable importance')
  )

  ggplotly(p, tooltip = sw)


### fall data only ####
# create empty list to hold dredges
dredges <- list()
# for each aphid taxon, make a global model and dredge through all predators
for (i in 1:length(short_aphlist)) {

  taxon <- short_aphlist[[i]]
  # note: mean_density_wide already has log-transformed arthropod data
  formula <-  paste0(taxon, '~ Arachnida + Coccinellidae + ',
                     'Ichneumonoidea + Nabis + Geocoris + ',
                     'Anthocoridae + (1|Site) + (1|Field)')
  mGlobal <- mean_density_wide %>%
    filter(Season=='Fall') %$%
    lmer(formula(formula),
         na.action = 'na.fail',
         REML = FALSE)

  dredges[[i]] <- dredge(mGlobal) # dredge and store output in list

}
names(dredges) <- short_aphlist

# check data sources
importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(ExpVars = names,
           VarWeights = sw,
           .keep = 'none') %>%
    arrange(ExpVars)
})

names(importance_tab) <- short_aphlist
# make tables of model selection results
FallTabs <- lapply(dredges, slice_head, n = 5)
names(FallTabs) <- short_aphlist
# show tables
# View(FallTabs)

# extract and examine the top Acyrthosiphon model.best.acy.mod <- get.models(FallTabs[[1]], subset = 1)[[1]]
summary(best.acy.mod)
# must remake to plot effects.
best.acy.mod.fa <- lmer(Acyrthosiphon ~ Anthocoridae + Ichneumonoidea +
                          (1|Site) + (1|Field),
                        data = mean_density_wide %>% filter(Season == 'Fall'),
                        REML = FALSE, na.action = 'na.omit')
summary(best.acy.mod.fa)

## better plot
# png('lattice.png', height = 6, width = 6, units = 'in', res = 300)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=2),
                     par.main.text = list(col = "blue", cex=0.5)))
plot(effect('Anthocoridae',best.acy.mod.fa,
            residuals = T),
     partial.residuals = list(smooth=F),
     axes = list(x = list(Anthocoridae = list(lab = 'Log(Anthocorid density)')),
                 y = list(lab = 'Log(Aphid density)')),
     main = NULL,
     lattice=list(key.args=list(axis.text = list(cex=10))))
plot(effect('Ichneumonoidea',best.acy.mod.fa,
            residuals = T),
     partial.residuals = list(smooth=F),
     axes = list(x = list(Ichneumonoidea = list(lab = 'Log(Ichneumonid density)')),
                 y = list(lab = 'Log(Aphid density)')),
     main = NULL,
     lattice=list(key.args=list(axis.text = list(cex=10))))
# dev.off()

# extract and examine the top Aphis model.
best.aphis.mod <- get.models(FallTabs[[2]], subset = 1)[[1]]
summary(best.aphis.mod)
# Must remake to plot effects.
best.aphis.mod <- lmer(Aphis ~ Arachnida + Ichneumonoidea + (1|Site) + (1|Field),
                       data = mean_density_wide %>% filter(Season == 'Fall'),
                       REML = FALSE)
# make plot
# png('Fall_aphis_effect.jpg',
#     width = 7,
#     height = 5,
#     units = 'in',
#     res = 300)
plot(allEffects(best.aphis.mod, residuals = TRUE),
     main = 'Aphis, Fall',
     id = list(n = 36, labels = mean_density_wide %>%
                 filter(Season == 'Fall') %>%
                 unite('id', Site, Field, Plot, sep = '.') %>%
                 pull(id)))
# dev.off()

# extract and examine the top Therioaphis model.
best.therio.mod <- get.models(FallTabs[[3]], subset = 1)[[1]]
summary(best.therio.mod)
# Must remake to plot effects.
best.therio.mod <- lmer(Therioaphis ~ Arachnida + (1|Site) + (1|Field),
                        data = mean_density_wide %>% filter(Season == 'Fall'),
                        REML = FALSE)
# make plot
# png('Fall_therio_effect.jpg',
#     width = 7,
#     height = 5,
#     units = 'in',
#     res = 300)
plot(allEffects(best.therio.mod, residuals = TRUE),
     main = 'Therioaphis, Fall',
     id = list(n = 36, labels = mean_density_wide %>%
                 filter(Season == 'Fall') %>%
                 unite('id', Site, Field, Plot, sep = '.') %>%
                 pull(id)))
# dev.off()

# for each aphid taxon, make variable importance heatmap for all aphids
aphPredVIFall <- list()
for (i in 1:length(short_aphlist)){

  p <- bind_rows(importance_tab, .id = 'Taxon') %>%
    mutate(VarWeights = as.numeric(VarWeights),
           ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
           Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
    filter(Taxon == short_aphlist[[i]]) %>%
    ggplot(aes(x = reorder(ExpVars, desc(VarWeights)),
               y = VarWeights,
               fill = VarWeights)) +
    geom_col() +
    scale_fill_gradient('',
                        low="blue",
                        high="red",
                        breaks = c(0.3,0.9),
                        labels = c('low','high')) +
    labs(title = short_aphlist[[i]], x = '', y = 'Variable importance') +
    theme(legend.position = c(0.95,0.85))

  aphPredVIFall[[i]] <- p

}
# print all plots
for (i in 1:length(aphPredVIFall)) {

  print(aphPredVIFall[[i]])

}
# make a larger heatmap that includes all aphid taxa in one figure

# combine all importance tables into one
aphVIFall <- rbindlist(importance_tab, idcol = TRUE) %>% as_tibble()
# constuct heatmap
PLOT_aphVIFall <- ggplot(data = aphVIFall,
                     aes(x = .id, y = ExpVars, fill = VarWeights)) +
  geom_tile() +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45))+
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Aphid',
       y = 'Predator',
       title = 'Variable importance, aphids ~ predators') +
  coord_flip()
# show it
PLOT_aphVIFall








# clean environment
rm(dredges, importance_tab, mGlobal, springTabs, FallTabs)

# Predators ~ ####
### Prepare data ####
# already done. see:
landCoverTabs

#### functions ####
# extract and examine models by rank
buildBestLandcoverMod <- function(mod_table, rank, season = 'unknown') {

  best.mod <- get.models(mod_table, subset = rank)[[1]]
  cat(green('trying rank', rank, '\n'))
  while (ncol(best.mod@frame) <= 2) {

    cat(yellow('rank', rank, 'model is null\n'))

    rank = rank + 1
    cat(green('trying rank', rank, '\n'))
    best.mod <- get.models(mod_table, subset = rank)[[1]]
  }
  getPred <- best.mod@call$formula[[2]]
  cat(green('Summary for model',  'rank =', rank, '\n'))
  # plot(allEffects(best.mod, residuals = TRUE),
  #      main = paste(paste(getPred, collapse = ''),
  #                   'landcover model rank',
  #                   as.character(rank),
  #                   'Season:', season),
  #      id = list(n = length(get(best.mod@call$data)$id),
  #                labels = get(best.mod@call$data)$id))
  summary(best.mod)
}

# plot the variable importance
plotVarImportance <- function (mod_table, season = 'unknown season'){

  # get mod tab name
  dataset <- (names(mod_table))
  print(dataset)
  mod_table <- mod_table[[1]]
  # identify the response var
  getPred <- get.models(mod_table, subset = 1)[[1]]@call$formula[[2]]


  importance_tab <- sw(mod_table) %>% #tibble(names = names(.))
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    arrange(names) %>%
    separate(names, c('class', 'distWeight'), sep = "_") %>%
    mutate(distWeight = as_factor(recode(distWeight,
                                         `sig1` = 'Very aggressive',
                                         `sig2` = 'Aggressive',
                                         `sig3` = 'Moderately aggressive',
                                         `sig4` = 'Moderate',
                                         `sig5` = 'Slight',
                                         `sig6` = 'Minimal',
                                         `const` = 'Constant',
                                         `no` = 'None'))) %>%
    mutate(distWeight = fct_relevel(distWeight, 'Constant', 'None', after = Inf),
           class = recode(class,
                          ag = 'Agricultural',
                          alfalfa = 'Alfalfa',
                          dirt = 'Bare soil +\n dirt road',
                          impermeable = 'Impermeable',
                          naturalArid = 'Natural',
                          water = 'Surface\nwater',
                          weedy = 'Weedy',
                          wet = 'Riparian'))

  p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
    geom_tile() +
    theme(
      axis.text.x=element_text(angle = 45, hjust = 0),
      axis.text.y = element_text(angle = 45))+
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         title = paste0('log(', getPred[2], ')',
                       ' Variable importance, ',
                       season)
         )

  pp <- ggplotly(p, tooltip = 'sw')

  group_importance <- importance_tab %>%
    group_by(distWeight) %>%
    summarize(weight = sum(sw))

  q <- ggplot(data = group_importance,
              aes(x = '', y = distWeight, fill = weight)) +
    geom_tile() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         fill = '')

  qq <- ggplotly(q, tooltip = 'weight')
  subplot(pp, qq, widths = c(7/8, 1/8)) %>%
    layout(title = dataset)

}

# Add models with margin data to model selection table.
makeGlobalLandcoverModTabSpring <- function(mod_table, rank = 1, m.max = 5){
  # mod_table = arachnida_mod_table_sp

 mod <- get.models(mod_table, subset = rank)[[1]]
  # Global model is rank-deficient. We will have to 'trick' dredge.

  # mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
  #                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
  #                     (1|Site),
  #                  data = fD_spring_sub) # This doesn't work.
  response <- names(mod@frame) %>% head(1)
  landCoverPredictors <- names(mod@frame) %>% head(-1) %>% tail(-1)
  tryCatch({

    distWeight <- str_match(landCoverPredictors[[1]], '_[:alnum:]+')

  }, error = function(e){
    print(e)
    distWeight <<- '_no'
  })
  predictors <- paste0(varList,
                       distWeight)
  vars.all <- c(predictors,
                'shan',
                'rich',
                'total_cover')
  # Write formula for full global model.
  form <- formula(paste0(response, '~',
                         paste0(vars.all, collapse='+'),
                         '+(1|Site)'))

  # Prepare formula with reduced number of terms.
  form.red <- formula(paste0(response, ' ~ (1|Site)')) #error was here
  # Fit reduced model.
  fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE,
                   na.action = 'na.fail')
  # Replace reduced model formula with full global model formula.
  attr(fmod.red@frame, "formula") <- form
  # Check formula attribute of reduced model.
  formula(fmod.red) # Looks good.

  # Run dredge() with m.max parameter to avoid convergence failures.
  ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = m.max), eval))
  ms1
}

makeGlobalLandcoverModTabFall <- function(mod_table, rank = 1, m.max = m.max){
  # mod_table = arachnida_mod_table_sp

  mod <- get.models(mod_table, subset = rank)[[1]]
  # Global model is rank-deficient. We will have to 'trick' dredge.

  # mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
  #                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
  #                     (1|Site),
  #                  data = fD_spring_sub) # This doesn't work.
  response <- names(mod@frame) %>% head(1)
  landCoverPredictors <- names(mod@frame) %>% head(-1) %>% tail(-1)
  tryCatch({

    distWeight <- str_match(landCoverPredictors[[1]], '_[:alnum:]+')

  }, error = function(e){
    print(e)
    distWeight <<- '_no'
  })
  predictors <- paste0(varList,
                       distWeight)
  vars.all <- c(predictors,
                'shan',
                'rich',
                'total_cover')
  # Write formula for full global model.
  form <- formula(paste0(response, '~',
                         paste0(vars.all, collapse='+'),
                         '+(1|Site)'))

  # Prepare formula with reduced number of terms.
  form.red <- formula(paste0(response, ' ~ (1|Site)')) #error was here
  # Fit reduced model.
  fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE,
                   na.action = 'na.fail')
  # Replace reduced model formula with full global model formula.
  attr(fmod.red@frame, "formula") <- form
  # Check formula attribute of reduced model.
  formula(fmod.red) # Looks good.

  # Run dredge() with m.max parameter to avoid convergence failures.
  ms1 <- model.sel(lapply(dredge(fmod.red,
                                 evaluate = FALSE,
                                 m.max = m.max), eval))
  ms1
}


# Plot the variable importance.
plotGlobalVarImportance <- function(global_tab){
  # global_tab <- araG1
  importance_tab <- sw(global_tab) %>%
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    separate(names, c('names', 'weights'), sep = '_') %>%
    mutate(names = fct_recode(names,
                              'Agricultural'='ag',
                              'Alfalfa'='alfalfa',
                              'Bare soil +\n dirt road'='dirt',
                              'Impermeable'='impermeable',
                              'Natural'='naturalArid',
                              'Surface\nwater'='water',
                              'Weedy'='weedy',
                              'Riparian'='wet',
                              'Shannon D' = 'shan',
                              'Total\ncover'='total',
                              'Richness'='rich')) %>%
    rename(varWeight = sw)

  # function(rx) str_view_all("bacad", rx)


  p <- ggplot(data = importance_tab,
              aes(x = reorder(names,-varWeight), y = varWeight)) +
    geom_col() +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
    labs(x = 'Variable', y = 'Importance')

  ggplotly(p, tooltip = 'sw')
  p
}

# Build final model based on variable selection exercise.
buildFinalLandcoverMod <- function(no_margin_tab, global_tab, rank = 1){

  gMod <- get.models(global_tab, subset = rank)[[1]]


  while (ncol(gMod@frame) <= 2) {


    cat(yellow('global mod rank', rank, 'model is null\n'))

    rank = rank + 1
    cat(green('trying rank', rank, '\n'))
    gMod <- get.models(global_tab, subset = rank)[[1]]
  }

  globalPredictors <- names(gMod@frame) %>% head(-1) %>% tail(-1)

  rank = 1
  nmMod <- get.models(no_margin_tab, subset = rank)[[1]]

  while (ncol(nmMod@frame) <= 2) {

    cat(yellow('no-margin mod rank', rank, 'model is null\n'))

    rank = rank + 1
    cat(green('trying rank', rank, '\n'))
    nmMod <- get.models(no_margin_tab, subset = rank)[[1]]

  }


  noMarginPredictors <- names(nmMod@frame) %>% head(-1) %>% tail(-1)

  response <- names(gMod@frame) %>% head(1)

  if (TRUE %in% (globalPredictors %in% c('rich', 'shan', 'total_cover'))){

    cat(yellow('Margin data included in model;\nDropping Yerington data'))
    fin.mod <- lmer(as.formula(paste0(response, '~', globalPredictors,'+ (1|Site)')),
                    data = fD_spring_sub, REML = FALSE)
  } else {
    cat(yellow('No margin data necessary;\nusing all sites'))
    fin.mod <- nmMod
  }
  getPred <- fin.mod@call$formula[[2]]
  cat(green('Summary for final model',  'rank =', '??', '\n'))
  plot(allEffects(fin.mod, residuals = TRUE),
       main = paste(paste(getPred, collapse = ''),
                    'landcover model rank',
                    as.character(rank),
                    'Season:', '???'),
       id = list(n = length(get(fin.mod@call$data)$id),
                 labels = get(fin.mod@call$data)$id))
  summary(fin.mod)
  fin.mod

}


### Correlation of landcover factors ####
# list of distance weights (already made)
weightList
# new names for distance weights
distList <- c('Very Aggressive',
              'Aggressive',
              'Moderately aggressive',
              'Moderate',
              'Slight',
              'Constant',
              'None')
# make list to hold ....
dVarList <- list()
# for each dataset ...
for (i in 1:length(landCoverTabs)){

  # for each distance weight...
  tabList <- list()
  for (j in 1:length(distList)) {

    # build correlation tables
    if (i %% 2 == 0){ # if landcover tab is 8 class:
    tabList[[j]] <- landCoverTabs[[i]] %>%
      select(ends_with(weightList[[j]])) %>%
      rename_with(~ str_remove(.x, "_.+")) %>%
      rename('Agricultural' = ag,
             'Alfalfa' = alfalfa,
             'Bare soil +\n dirt road' = dirt,
             'Impermeable' = impermeable,
             'Natural' = naturalArid,
             'Surface\nwater' = water,
             'Weedy' = weedy,
             'Riparian' = wet) %>%
      relocate(sort(peek_vars())) %>%
      correlate %>%
      shave
    } else { # if landcover tab is 7 class:
      tabList[[j]] <- landCoverTabs[[i]] %>%
        select(ends_with(weightList[[j]])) %>%
        rename_with(~ str_remove(.x, "_.+")) %>%
        rename('Agricultural' = ag,
               'Alfalfa' = alfalfa,
               'Bare soil +\n dirt road' = dirt,
               'Impermeable' = impermeable,
               'Natural' = naturalArid,
               'Surface\nwater' = water,
               'WeedyRiparian' = weedyWet) %>%
        relocate(sort(peek_vars())) %>%
        correlate %>%
        shave

    }

  }
  names(tabList) <- distList
  dVarList[[i]] <- tabList
}
names(dVarList) <- names(landCoverTabs)
# # check
# dVarList

# build correlation plots
# for each data set...
PLOTS_corr <- list()
for (i in 1:length(dVarList)) {

  # for each distance weight...
  tempList <- list()
  for (j in 1:length(distList)) {
    p <- rplot(dVarList[[i]][[j]],
                print_cor = TRUE,
                legend = FALSE,
                .order = 'alphabet') +
            theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
            labs(title = paste(names(dVarList[i]),distList[[j]]))
    tempList[[j]] <- p


  }
  names(tempList)<-distList
  PLOTS_corr[[i]] <- tempList
}
names(PLOTS_corr) <- names(landCoverTabs)

### Flooded vs sprinklers ####
# make df of watering methods
wateringMethod <- c(
  'Flooding',
  'Flooding',
  'Flooding',
  'Flooding',
  'Flooding',
  'Flooding',
  'Sprinklers',
  'Sprinklers',
  'Sprinklers',
  'Sprinklers',
  'Flooding',
  'Flooding'
)
idList <- landCoverTabs[[1]] %>% filter(Season == 'Spring') %>%
  pull(id)
wM <- tibble(wateringMethod, idList)
wM
# join this new df to the arthropod data
withWM <- mean_density_wide %>%
  mutate(idList = paste0(Site, '0', Field)) %>%
  left_join(wM, by = c('idList' = 'idList'))

# model
floodMod <- lmer(Acyrthosiphon ~ wateringMethod + (1|Site),
                 data = withWM)
floodMod
summary(floodMod)

plot(allEffects(floodMod))

# some plots
p <- ggplot(withWM,
       aes(x = wateringMethod,
           y = Acyrthosiphon)) +
  # geom_boxplot() +
  # geom_jitter() +
  # geom_violin()
  geom_jitter(aes(color = Site), width = 0.06) +
  # geom_smooth(method = 'lm') +
  stat_summary(fun.data = 'mean_cl_boot') +
  labs(title = 'mean_cl_boot',
       subtitle = paste0('basic nonparametric bootstrap for ',
         'obtaining confidence limits')) +
  scale_color_brewer(palette = 'Set1')

p
q=ggplot(withWM,
       aes(x = wateringMethod,
           y = Acyrthosiphon)) +
  # geom_boxplot() +
  # geom_jitter() +
  geom_violin() +
  geom_jitter(aes(color = Site), width = 0.06) +
  # geom_smooth(method = 'lm') +
  stat_summary(fun.data = 'mean_cl_normal') +
  labs(title = 'mean_cl_normal',
       subtitle = paste0('lower and upper Gaussian confidence limits ',
         'based on the t-distribution')) +
  scale_color_brewer(palette = 'Set1')
q

r=ggplot(withWM,
       aes(x = wateringMethod,
           y = Acyrthosiphon)) +
  geom_boxplot() +
  # geom_jitter() +
  # geom_violin()
  geom_jitter(aes(color = Site), width = 0.06) +
  # geom_smooth(method = 'lm') +
  # stat_summary(fun.data = 'mean_cl_normal') +
  labs(title = 'boxplot',
       subtitle = 'a normal boxplot') +
  scale_color_brewer(palette = 'Set1')
r
grid.arrange(p,q,r, nrow =1)

# t tests
flooded <- withWM %>%
  filter(wateringMethod=='Flooding') %>%
  select(Acyrthosiphon)
sprinklers <- withWM %>%
  filter(wateringMethod=='Sprinklers') %>%
  select(Acyrthosiphon)

t.test(flooded, sprinklers, var.equal = TRUE)
t.test(flooded, sprinklers)

final_fig <- ggplot(withWM,
                    aes(x = wateringMethod,
                        y = Coccinellidae)) +
  geom_jitter(aes(color = Site), width = 0.1) +
  stat_summary(fun.data = 'mean_cl_boot', geom="errorbar", width = 0.3)+
  labs(
       y = 'log(Ladybug density)',
       x = 'Irrigation method')+
  # theme_classic(base_size = 20) +
  scale_color_brewer(palette = 'Set1')
final_fig
# ggsave('watering.png', width = 5.5, height = 6.5)



# ### sum of nat ####
# # note: this code is for a demo figure in the ESA presentation
# sums <- field_data %>%
#   filter(Season == 'Spring') %>% # must not use scaled data
#   mutate(sum_no = rowSums(across(ends_with('no')))) %>%
#   mutate(sum_sig4 = rowSums(across(ends_with('sig4')))) %>%
#   mutate(sum_sig2 = rowSums(across(ends_with('sig2')))) %>%
#   select(id, ends_with('no'), ends_with('sig4'), ends_with('sig2')) %>%
#   select(id, starts_with('naturalArid'), starts_with('sum')) %>%
#   mutate(ratio_no = (naturalArid_no/sum_no)*100,
#          ratio_sig4 = (naturalArid_sig4/sum_sig4)*100,
#          ratio_sig2 = (naturalArid_sig2/sum_sig2)*100) %>%
#   select(id, starts_with('ratio')) %>%
#   filter(id == 'Minden02')
#
#
# sums
# ratios<-sums %>% select(-id)
# percentile <- ratios$ratio_no
# ggplot() +
#   geom_col(aes("", 100)) +
#   geom_col(aes("", percentile), fill = "#CF3FFF") +
#   coord_flip() +
#   theme_void()
# ggsave('naNo.png', width = 3, height = 0.5, units = 'in')
# percentile <- ratios$ratio_sig4
# ggplot() +
#   geom_col(aes("", 100)) +
#   geom_col(aes("", percentile), fill = "#CF3FFF") +
#   coord_flip() +
#   theme_void()
# ggsave('sig4No.png', width = 3, height = 0.5, units = 'in')
# percentile <- ratios$ratio_sig2
# ggplot() +
#   geom_col(aes("", 100)) +
#   geom_col(aes("", percentile), fill = "#CF3FFF") +
#   coord_flip() +
#   theme_void()
# ggsave('sig2No.png', width = 3, height = 0.5, units = 'in')


#### General figs ####
# make list of all rds files
# RDS filename format:
# Geocoris_fixed8_sub_fall
# taxaList_dataset_yeringtonIncluded?_season

# make list of RDS filenames
taxaList <- c('Arachnida',
           'Coccinellidae',
           'Geocoris',
           'Ichneumonoidea',
           'Acyrthosiphon')
datasetList <- c('regular7', 'regular8', 'fixed7', 'fixed8')
fS <- c('full', 'sub')
ssn <- c('spring', 'fall')
rdsList <- expand.grid(taxaList, datasetList, fS, ssn) %>%
  tibble %>%
  unite(rdsName)
# rdsList
# str(rdsList$rdsName)

# read in tabs
# note: this uses CRAZY ram
tabList <- c()
for (i in 1:length(rdsList$rdsName)) {

  tabList[[i]] <- readRDS(paste0('modTabs/',rdsList$rdsName[[i]]))

}
names(tabList) <- rdsList$rdsName

# # example grep filtering of modTabs
# tabList[grep("spring", names(tabList))]

# old
### spring data only

# # Subset data for modeling with landcover classes.
# fD_spring <- landCoverTabs[[1]] %>%
#   filter(Season == 'Spring') %>%
#   mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.
#
# # Subset data for modeling with landcover classes and margin data.
# fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')

# generate vi heatmaps ####
# filter to focus on relevant taxa
taxaTabList <- tabList[grep(paste(taxaList,collapse='|'),
                                  names(tabList))]
# split tablist by season
springTabList <- taxaTabList[grep("spring", names(taxaTabList))]
fallTabList <- taxaTabList[grep("fall", names(taxaTabList))]
ssnList <- list(springTabList, fallTabList)
names(ssnList) <- c('spring', 'fall')

# make heatmaps
PLOTS_heatmaps <- list()
for (i in 1:length(ssnList)) {

  tempList <- ssnList[[i]]
  viList <- c()
  for (j in 1:length(tempList)) {

    viList[[j]] <- plotVarImportance(tempList[j], names(ssnList[1]))

  }
  names(viList) <- names(tempList)
  PLOTS_heatmaps[[i]] <- viList
}
names(PLOTS_heatmaps) <- names(ssnList)

# 4 relevant plots for each taxon, using only "full" data
# group by taxon
PLOTS_4byHeatmaps <- list()
for (i in 1:length(PLOTS_heatmaps)) {

  tempList <- PLOTS_heatmaps[[i]]
  PLOTS_VI <- list()
  for (j in 1:length(taxaList)) {

    taxon <- taxaList[[j]]
    all <- tempList[grep(taxon, names(tempList))]
    onlyFull <- all[grep("full", names(all))]
    PLOTS_VI[[j]] <- onlyFull

  }
  names(PLOTS_VI) <- taxaList
  PLOTS_4byHeatmaps[[i]] <- PLOTS_VI
}
names(PLOTS_4byHeatmaps) <- names(ssnList)

# ugly plots
PLOTS_by4 <- list()
for (i in 1:length(ssnList)) {

  tempList <- PLOTS_4byHeatmaps[[i]]
  PLOTS_VI4panel <- list()
  for (j in 1:length(taxaList)) {

    PLOTS_VI4panel[[j]] <- subplot(tempList[[j]], nrows = 2) %>%
      layout(title = paste0(names(tempList[j]),
                            names(ssnList[i]),
                            ' <-7', # 7-class data on the left
                            # reg (not fixed) alfalfa class on the top
                            ' ^reg. classification'))

  }
  names(PLOTS_VI4panel) <- taxaList
  PLOTS_by4[[i]] <- PLOTS_VI4panel
}
names(PLOTS_by4) <- names(ssnList)


# get best models ####
# collect best models from each modTable
bestModsList <- list()
# for each season:
for (i in 1:length(ssnList)) {
  # grab a set of modTabs
  tempList <- ssnList[[i]]
  # to hold best models for each taxon
  bestMods <- list()
  # for each taxon
  for (j in 1:length(taxaList)) {
    # get the taxon
    taxon <- taxaList[[j]]
    # get the modTabs for that taxon
    all <- tempList[grep(taxon, names(tempList))]
    # list to hold all versions of models for a given taxon and season
    modList <- list()
    for (k in 1:length(all)) {
      # get the best non-null model
      modList[[k]] <- buildBestLandcoverMod(all[[k]], 1)

    }
    names(modList) <- names(all)
    bestMods[[j]] <- modList

  }
  names(bestMods) <- taxaList
  bestModsList[[i]] <- bestMods
}
names(bestModsList) <- names(ssnList)


# view head of tabs
heads <- c()
for (i in 1:length(tabList)) {

  heads[[i]] <- head(tabList[[i]])

}
names(heads) <- names(tabList)

## ladybug spring summaries
temp <- bestModsList$spring$Coccinellidae
for (i in 1:4){
  cat(red(names(temp)[i]))
  print(summary(temp[[i]]))
}

# plot top mod effects ####
# looping
# won't work bc of the way mods were built

# plot top mod effects (one at a time, manually)
# spring ladybugs first
temp <- ssnList$spring
# reg 7
reg7full <- temp$Coccinellidae_regular7_full_spring
plotMod <- get.models(reg7full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7 %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Coccinellidae ~ alfalfa_sig1 + dirt_sig1 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('reg7'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# reg 8
reg8full <- temp$Coccinellidae_regular8_full_spring
plotMod <- get.models(reg8full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8 %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Coccinellidae ~ dirt_sig1 + weedy_sig1 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('reg8'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# fix 7
fix7full <- temp$Coccinellidae_fixed7_full_spring
plotMod <- get.models(fix7full, 8)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7Fixed %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Coccinellidae ~ alfalfa_no + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('fix7'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))
# fix 8
fix8full <- temp$Coccinellidae_fixed8_full_spring
plotMod <- get.models(fix8full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8Fixed %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Coccinellidae ~ weedy_sig2 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('fix8'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# acyrthosiphon
reg7full <- temp$Acyrthosiphon_regular7_full_spring
plotMod <- get.models(reg7full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7 %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Acyrthosiphon ~ water_sig5 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('reg7'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# reg 8
reg8full <- temp$Acyrthosiphon_regular8_full_spring
plotMod <- get.models(reg8full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8 %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Acyrthosiphon ~ dirt_sig1 + weedy_sig1 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('reg8'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# fix 7
fix7full <- temp$Acyrthosiphon_fixed7_full_spring
plotMod <- get.models(fix7full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7Fixed %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Acyrthosiphon ~ water_sig5 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('fix7'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

# fix 8
fix8full <- temp$Acyrthosiphon_fixed8_full_spring
plotMod <- get.models(fix8full, 1)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8Fixed %>% filter(Season == 'Spring')
# remake with local data for effects plots
localPlotMod <- lmer(Acyrthosiphon ~ water_sig5 + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
cat(red('fix8'))
summary(plotMod)
summary(localPlotMod)
plot(allEffects(localPlotMod, residuals = TRUE))

## acyrthosiphon full models for spring ####
# irrigation data
wM
# join to datasets and filter to spring only
lcTabsWater <- list()
for (i in 1:length(landCoverTabs)){

  lcTabsWater[[i]] <- left_join(landCoverTabs[[i]],
                                wM,
                                by=(c('id'='idList'))) %>%
    filter(Season == 'Spring')

}
names(lcTabsWater) <- names(landCoverTabs)

# new version of modTab builder specifically for Acyrthosiphon


acySpTabs <- list()
for (j in 1:length(lcTabsWater)){

  data <- lcTabsWater[[j]]
  taxon='Acyrthosiphon'
  m.max=3

  distList <- c('_no ',
                '_const ',
                '_sig1 ',
                '_sig2 ',
                '_sig3 ',
                '_sig4 ',
                '_sig5 ')

  varList <- data %>%
    select(contains('_')) %>%
    select(-total_cover) %>% # must remove this or it breaks
    names() %>% str_extract('[:alpha:]+') %>%
    unique()



  cand_mod_tabs <- list()

  if (taxon == 'empty' | !is_tibble(data)) {
    stop(red("Please specify taxon and data \n"), call. = FALSE)
  }

  # print inputs
  cat(yellow('Taxon:'),
      green(taxon),
      yellow('Data:'),
      green(deparse(substitute(data))),
      '\n')

  # get dfname
  dfname <- as.name(deparse(substitute(data)))
  # Build models
  for (i in 1:length(distList)) {
    # i=2

    # incase full model fails to fit, try 'tricking' dredge
    message(blue('fitting', distList[[i]],'models'))
    # fit reduced model
    fmod.red <- lmer(as.formula(
      paste0(taxon, ' ~ (1|Site)')),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')
    # define full model formula
    form <-formula(
      paste0(taxon, ' ~ ',
             paste0(varList, distList[[i]], '+ ', collapse = ''),
             'Coccinellidae + wateringMethod + (1|Site)'  # add predators and water
      ))
    # Replace reduced model formula with full global model formula.
    attr(fmod.red@frame, "formula") <- form
    # Run dredge() with m.max parameter to avoid convergence failures.
    fmod.red@call$data <- dfname

    cand_mod_tabs[[i]] <-  # superassign?
      model.sel(lapply(
        dredge(fmod.red, m.lim = c(NA, m.max), trace = 2, evaluate = FALSE),
        eval))

    message(blue(nrow(cand_mod_tabs[[i]]), 'models fit\n'))

  }
  # Rbind the elements of the list together. This forces recalculation of AICc
  mod_table <- rbind(cand_mod_tabs[[1]],
                     cand_mod_tabs[[2]],
                     cand_mod_tabs[[3]],
                     cand_mod_tabs[[4]],
                     cand_mod_tabs[[5]],
                     cand_mod_tabs[[6]],
                     cand_mod_tabs[[7]])

  acySpTabs[[j]] <- mod_table
}

names(acySpTabs) <- names(lcTabsWater)
print(acySpTabs[[1]])
# acy spring specific heatmap func
plotVarImportance2 <- function (mod_table,
                                season = 'unknown season'){
  # mod_table = acySpTabs[1]
  season = 'Spring'
  # get mod tab name
  dataset <- (names(mod_table))
  print(dataset)
  mod_table <- mod_table[[1]]



  # identify the response var
  getPred <- 'Acyrthosiphon'


  importance_tab <- sw(mod_table) %>% #tibble(names = names(.))
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    arrange(names) %>%
    separate(names, c('class', 'distWeight'), sep = "_") %>%
    mutate(distWeight = as_factor(recode(distWeight,
                                         `sig1` = 'Very aggressive',
                                         `sig2` = 'Aggressive',
                                         `sig3` = 'Moderately aggressive',
                                         `sig4` = 'Moderate',
                                         `sig5` = 'Slight',
                                         `sig6` = 'Minimal', # wtf did this come from
                                         `const` = 'Constant',
                                         `no` = 'None'))) %>%
    mutate(distWeight = fct_relevel(distWeight, 'Constant', 'None', after = Inf),
           class = recode(class,
                          ag = 'Agricultural',
                          alfalfa = 'Alfalfa',
                          dirt = 'Bare soil +\n dirt road',
                          impermeable = 'Impermeable',
                          naturalArid = 'Natural',
                          water = 'Surface\nwater',
                          weedy = 'Weedy',
                          wet = 'Riparian',
                          weedyWet = 'Weedy +\n Riparian'))

  p <- ggplot(data = importance_tab,
              aes(x = fct_relevel(class, 'Coccinellidae', 'wateringMethod', after = Inf),
                  y = distWeight,
                  fill = sw)) +
    geom_tile() +
    theme(
      axis.text.x=element_text(angle = 45, hjust = 0),
      axis.text.y = element_text(angle = 45))+
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         title = paste0('log(', getPred, ')',
                        ' Variable importance, ',
                        season)
    )
  p
  pp <- ggplotly(p, tooltip = 'sw')

  group_importance <- importance_tab %>%
    group_by(distWeight) %>%
    summarize(weight = sum(sw))

  q <- ggplot(data = group_importance,
              aes(x = '', y = distWeight, fill = weight)) +
    geom_tile() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         fill = '')

  qq <- ggplotly(q, tooltip = 'weight')
  subplot(pp, qq, widths = c(7/8, 1/8)) %>%
    layout(title = dataset)

}

acyHeatMaps <- list()
for (i in 1:4){
acyHeatMaps[[i]]<-plotVarImportance2(acySpTabs[i], 'Spring')
print(acyHeatMaps[[i]])
}





# check margin effects ####
# ladybugs
# reg 7
reg7sub <- temp$Coccinellidae_regular7_sub_spring
plotMod <- get.models(reg7sub, 8)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7 %>% filter(Season == 'Spring',
                                                Site != 'Yerington')
# remake with local data for effects plots
# using data from best "full" mod
localPlotMod <- lmer(Coccinellidae ~ dirt_sig1 + alfalfa_sig1 +
                       shan + rich + total_cover + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
margTab <- dredge(localPlotMod, m.max = 3)

importance_tab <- sw(margTab) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names)

ggplot(data = importance_tab, aes(x = reorder(names, sw), y = sw, fill = sw)) +
  geom_col() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Factor',
       y = 'Weight')

# reg 8
reg8sub <- temp$Coccinellidae_regular8_sub_spring
plotMod <- get.models(reg8sub, 8)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8 %>% filter(Season == 'Spring',
                                                Site != 'Yerington')
# remake with local data for effects plots
# using data from best "full" mod
localPlotMod <- lmer(Coccinellidae ~ dirt_sig1 + weedy_sig1 +
                       shan + rich + total_cover + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
margTab <- dredge(localPlotMod, m.max = 3)

importance_tab <- sw(margTab) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names)

ggplot(data = importance_tab, aes(x = reorder(names, sw), y = sw, fill = sw)) +
  geom_col() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Factor',
       y = 'Weight')

# fix 7
fix7sub <- temp$Coccinellidae_fixed7_sub_spring
plotMod <- get.models(fix7sub, 8)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover7Fixed %>% filter(Season == 'Spring',
                                                Site != 'Yerington')
# remake with local data for effects plots
# using data from best "full" mod
localPlotMod <- lmer(Coccinellidae ~ alfalfa_no +
                       shan + rich + total_cover + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
margTab <- dredge(localPlotMod, m.max = 3)

importance_tab <- sw(margTab) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names)

ggplot(data = importance_tab, aes(x = reorder(names, sw), y = sw, fill = sw)) +
  geom_col() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Factor',
       y = 'Weight')

# fix 8
fix8sub <- temp$Coccinellidae_fixed8_sub_spring
plotMod <- get.models(fix8sub, 8)[[1]]
summary(plotMod)
localData <-landCoverTabs$landcover8Fixed %>% filter(Season == 'Spring',
                                                     Site != 'Yerington')
# remake with local data for effects plots
# using data from best "full" mod
localPlotMod <- lmer(Coccinellidae ~ weedy_sig2 +
                       shan + rich + total_cover + (1|Site),
                     data = localData,
                     REML = FALSE,
                     na.action = 'na.fail')
margTab <- dredge(localPlotMod, m.max = 3)

importance_tab <- sw(margTab) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names)

ggplot(data = importance_tab, aes(x = reorder(names, sw), y = sw, fill = sw)) +
  geom_col() +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Factor',
       y = 'Weight')

# quick ladybug plot
landCoverTabs$landcover8 %>% filter(Season == 'Spring') %>%
ggplot(aes(x = reorder(id, Acyrthosiphon), y = Acyrthosiphon, fill = Site)) +
  geom_col() +
  scale_fill_brewer(palette = 'Set1')





# plotMod <- get.models(acyG1, 4)[[1]]
summary(plotMod)
  nullMod <- lmer(log(Coccinellidae + 1)~(1|Site), data = fD_spring, REML = F)
  anova(plotMod, nullMod)
  png(paste0('fixeda',plotName),
      width = 6,
      height = 6,
      units = 'in',
      res = 300)
  plot(effect('wet_sig2',plotMod, residuals = TRUE),
       main = NULL,
       partial.residuals = list( smooth = F),
       xlab = 'Riparian cover,\nslight weighting',
       ylab = 'log(Aphid density)',
       id = list(n = length(get(plotMod@call$data)$id),
                 labels = get(plotMod@call$data)$id))
  dev.off()

  png(paste0('fixedb',plotName),
      width = 6,
      height = 6,
      units = 'in',
      res = 300)
  plot(effect('weedy_sig1',plotMod, residuals = TRUE),
       main = NULL,
       partial.residuals = list( smooth = F),
       xlab = 'Weedy cover,\nvery aggressive weighting',
       ylab = 'log(Ladybug density)',
       id = list(n = length(get(plotMod@call$data)$id),
                 labels = get(plotMod@call$data)$id))
  dev.off()



ggplot(fD_spring_sub, aes(reorder(id, shan, desc), shan)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = 'Field', y = 'Shannon diversity of field margin')
ggsave('shan.png', width = 12.5, height = 6, units = 'in')
# make global modsel tabs
# ara
araG1 <- makeGlobalLandcoverModTabSpring(tabList[[1]], 1, 3)
araG2 <- makeGlobalLandcoverModTabSpring(tabList[[1]], 2, 3)
plotGlobalVarImportance(araG1)
ggsave('araG1.png', width = 6.25, height = 6, units = 'in')
plotGlobalVarImportance(araG2)
ggsave('araG2.png', width = 6.25, height = 6, units = 'in')


ggplot(field_data, aes(x = reorder(id, wet_sig1), y = wet_sig1)) +
  geom_col() +
  labs(y = 'Riparian, very aggressive weighting',
       x = 'Field') +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave('riparian_sig1.png', width = 6, height = 4, units = 'in')
# cocc ####
coccG1 <- makeGlobalLandcoverModTabSpring(tabList[[2]], 1, 3)
coccG2 <- makeGlobalLandcoverModTabSpring(tabList[[2]], 5, 3)
plotGlobalVarImportance(coccG1)
ggsave('coccG1.png', width = 12.5, height = 6, units = 'in')
plotGlobalVarImportance(coccG2)
ggsave('coccG2.png', width = 6.25, height = 6, units = 'in')


c1mod<-lmer(log(Coccinellidae+1)~weedy_sig1+(1|Site),
            data = fD_spring, REML = FALSE)
c2mod<-lmer(log(Coccinellidae+1)~rich+(1|Site),
            data = fD_spring_sub, REML = FALSE)
c3mod<-lmer(log(Coccinellidae+1)~total_cover+(1|Site),
            data = fD_spring_sub, REML = FALSE)
cmods<-list(c1mod, c2mod, c3mod)
model.sel(cmods)

anova(c1mod, c2mod)
plot(allEffects(c1mod, residuals = T))

png(paste0('cmods.png'),
    width = 6.5,
    height = 6,
    units = 'in',
    res = 300)

plot(effect('weedy_sig1',c1mod, residuals = TRUE),
     main = NULL,
     ylab = 'log(Ladybug density)',
     xlab = 'Weedy cover,\n very aggressive weighting',
     id = list(n = length(get(c1mod@call$data)$id),
               labels = get(c1mod@call$data)$id))
dev.off()

png(paste0('cmods2.png'),
    width = 6.5,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('rich',c2mod, residuals = TRUE),
     main = NULL,
     xlab = 'Plant species richness in field margins',
     id = list(n = length(get(c2mod@call$data)$id),
               labels = get(c2mod@call$data)$id))
dev.off()
png(paste0('cmods3.png'),
    width = 6.5,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('total_cover',c3mod, residuals = TRUE),
     main = NULL,
     xlab = 'Plant cover in field margins',
     id = list(n = length(get(c2mod@call$data)$id),
               labels = get(c2mod@call$data)$id))
dev.off()
dev.off()
plot(allEffects(c2mod, residuals = T))
plot()

gplot(field_data, aes(x = reorder(id, weedy_sig4), y = weedy_sig4)) +
  geom_col() +
  labs(y = 'Weedy, moderate weighting',
       x = 'Field') +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave('weedy_sig4.png', width = 6, height = 4, units = 'in')




# geoc
geocG1 <- makeGlobalLandcoverModTabSpring(tabList[[3]], 8, 3)
geocG2 <- makeGlobalLandcoverModTabSpring(tabList[[3]], 9, 3)
plotGlobalVarImportance(geocG1)
ggsave('geocG1.png', width = 6.25, height = 6, units = 'in')
plotGlobalVarImportance(geocG2)
ggsave('geocG2.png', width = 6.25, height = 6, units = 'in')

#
View(geocG1) # top mod is null, mod2 is landcover only
View(geocG2) # top mod is null, mod2 is landcover only

geocGM1 <- get.models(geocG1, 2)[[1]]
png('geocGM1.png',
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(allEffects(geocGM1, residuals = TRUE),
     main = NULL,
     xlab = 'wet_sig1',
     id = list(n = length(get(geocGM1@call$data)$id),
               labels = get(geocGM1@call$data)$id))
dev.off()

geocGM2 <- get.models(geocG2, 4)[[1]]
png('geocGM2_2.png',
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(allEffects(geocGM2, residuals = TRUE),
     main = NULL,
     xlab = 'Richness in field margins',
     id = list(n = length(get(geocGM1@call$data)$id),
               labels = get(geocGM1@call$data)$id))
dev.off()


# ich
ichG1 <- makeGlobalLandcoverModTabSpring(tabList[[4]], 1, 5)
plotGlobalVarImportance(ichG1)
ggsave('ichG1.png', width = 12.5, height = 6, units = 'in')
ichG2 <- makeGlobalLandcoverModTabSpring(tabList[[4]], 8, 5)
plotGlobalVarImportance(ichG2)
ggsave('ichG2.png', width = 12.5, height = 6, units = 'in')

View(ichG1)
ichGM1 <- get.models(ichG1, 3)[[1]]
png('ichGM1_1.png',
    width = 6.25,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('impermeable_sig1',ichGM1, residuals = TRUE),
     main = NULL,
     xlab = 'Impermeable cover, very aggressive weighting',
     id = list(n = length(get(ichGM1@call$data)$id),
               labels = get(ichGM1@call$data)$id))
dev.off()
png('ichGM1_2.png',
    width = 6.25,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('total_cover',ichGM1, residuals = TRUE),
     main = NULL,
     xlab = 'Plant cover in field margins',
     id = list(n = length(get(ichGM1@call$data)$id),
               labels = get(ichGM1@call$data)$id))
dev.off()

View(ichG2)
ichGM2 <- get.models(ichG2, 3)[[1]]
png('ichGM2.png',
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('shan',ichGM2, residuals = TRUE),
     main = NULL,
     xlab = 'Shannon diversity in field margins',
     id = list(n = length(get(ichGM1@call$data)$id),
               labels = get(ichGM1@call$data)$id))
dev.off()

# acy
acyG1 <- makeGlobalLandcoverModTabSpring(tabList[[5]], 1, 3)
plotGlobalVarImportance(acyG1)
ggsave('acyG1.png', width = 6.25, height = 6, units = 'in')
acyG2 <- makeGlobalLandcoverModTabSpring(tabList[[5]], 3, 3)
plotGlobalVarImportance(acyG2)
ggsave('acyG2.png', width = 6.25, height = 6, units = 'in')

View(acyG1) #1 weedy5, 2 null, 3 rich, 4 shan
View(acyG2) #1 null, 2 rich, 3 impermeable3, 4 weedy3, 5 shan

acyGM1 <- get.models(acyG1, 3)[[1]]
acyGM2 <- get.models(acyG1, 4)[[1]]
# acyGM3 <- get.models(acyG2, 2)[[1]] # same as acyGM1
png('acyGM1.png',
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(allEffects(acyGM1, residuals = TRUE),
     main = NULL,
     xlab = 'Richness in field margins',
     id = list(n = length(get(ichGM1@call$data)$id),
               labels = get(ichGM1@call$data)$id))
dev.off()

png('acyGM2.png',
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(allEffects(acyGM2, residuals = TRUE),
     main = NULL,
     xlab = 'Shannon diversity in field margins',
     id = list(n = length(get(ichGM1@call$data)$id),
               labels = get(ichGM1@call$data)$id))
dev.off()




#### Arachnida ####
# Build global model selection table
arachnida_mod_table_sp <- readRDS('modTabs/Arachnida_spring_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = arachnida_mod_table_sp, rank = 1, 'Spring season')
# Plot var importance charts
plotVarImportance(mod_table = arachnida_mod_table_sp)
# make global mod table (single distweight)
arachnida_global_sp <- makeGlobalLandcoverModTab(arachnida_mod_table_sp, 1)
# plot var importance for global mod table
plotGlobalVarImportance(arachnida_global_sp)
# build best model
arachnida_fin_mod_sp <- buildFinalLandcoverMod(arachnida_mod_table_sp, arachnida_global_sp, 1) # Error
summary(arachnida_fin_mod_sp)
plot(allEffects(arachnida_fin_mod_sp, residuals = TRUE), main = 'Arachnida, final spring model')


#### Coccinellidae ####
# Build global model selection table
coccinellidae_mod_table_sp <- readRDS('modTabs/Coccinellidae_spring_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_sp, rank = 1,'Spring')
test <- lmer(log(Coccinellidae+1)~total_cover+(1|Site), data = fD_spring_sub, REML = FALSE)
summary(test)
plot(allEffects(test, residuals=T))

# make nice plot
plotMod <- get.models(coccinellidae_mod_table_sp, 1)[[1]]
png('spring_cocc_weedy_effect.jpg',
    width = 6,
    height = 5,
    units = 'in',
    res = 300)
plot(allEffects(plotMod, residuals = TRUE),
     main = NULL,
     xlab = 'Weighted proportion of weedy cover',
     id = list(n = length(get(plotMod@call$data)$id),
               labels = get(plotMod@call$data)$id))
dev.off()

# try naturalArid for comparison to sentinel
subCoccTab <- coccinellidae_mod_table_sp %>% filter(
  !is.na(naturalArid_const) |
  !is.na(naturalArid_no) |
  !is.na(naturalArid_sig1) |
  !is.na(naturalArid_sig2) |
  !is.na(naturalArid_sig3) |
  !is.na(naturalArid_sig4) |
  !is.na(naturalArid_sig5)

) %>% head(1)

plot2mod <- lmer(log(Coccinellidae + 1)~total_cover+(1|Site), data = fD_spring_sub, REML = TRUE)
png('spring_cocc_totcOnly_effect.jpg',
    width = 7,
    height = 5,
    units = 'in',
    res = 300)
plot(effect('total_cover', plot2mod, residuals = TRUE),
     main = NULL,
     xlab = 'total cover alone, spring',
     id = list(n = length(get(plot2mod@call$data)$id),
               labels = get(plot2mod@call$data)$id))
dev.off()

fD_spring %>% ggplot(aes(x = reorder(id, Coccinellidae), y = log(Coccinellidae + 1))) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Field')
ggsave(filename = 'coccFieldDensity.jpg', width = 7, height = 5, units = 'in')

fD_spring %>% ggplot(aes(x = reorder(id, Coccinellidae), y = weedy_sig4)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Field')
ggsave(filename = 'weedy4.jpg', width = 7, height = 5, units = 'in')



# Plot var importance charts
plotVarImportance(mod_table = coccinellidae_mod_table_sp, "Spring")
ggsave(filename = 'Cocc_var_importance.jpg', width = 7, height = 5, units = 'in')



# make global mod table (single distweight)
coccinellidae_global_sp <- makeGlobalLandcoverModTab(coccinellidae_mod_table_sp, 1)
# plot var importance for global mod table
plotGlobalVarImportance(coccinellidae_global_sp)
# build best model
coccinellidae_fin_mod_sp <- buildFinalLandcoverMod(coccinellidae_mod_table_sp, coccinellidae_global_sp, 1) # Error
summary(coccinellidae_fin_mod_sp)
plot(allEffects(coccinellidae_fin_mod_sp, residuals = TRUE), main = 'Coccinellidae, final spring model')

##### Compare old and new ####
# 2022-07-28T17:39:44Z this is old now bc of superdove data. Deleted.

#### Ichneumonidae ####
# Build global model selection table
ichneumonidae_mod_table_sp <- readRDS('modTabs/Ichneumonidae_8class') # not built yet
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = ichneumonidae_mod_table_sp, rank = 1,'Spring')
# Plot var importance charts
plotVarImportance(mod_table = ichneumonidae_mod_table_sp)
# make global mod table (single distweight)
ichneumonidae_global_sp <- makeGlobalLandcoverModTab(ichneumonidae_mod_table_sp, 1)
# plot var importance for global mod table
plotGlobalVarImportance(ichneumonidae_global_sp)
# build best model
ichneumonidae_fin_mod_sp <- buildFinalLandcoverMod(ichneumonidae_mod_table_sp,
                                         ichneumonidae_global_sp,
                                         2)
summary(ichneumonidae_fin_mod_sp)
plot(allEffects(ichneumonidae_fin_mod_sp, residuals = TRUE), main = 'Ichneumonidae, final spring model')
# hypothetical ich model with dirt_const and shan:
ich_test_mod_sp <- lmer(log(Ichneumonidae + 1) ~ dirt_const + shan + (1|Site),
                     data = fD_spring_sub,
                     REML = FALSE)
summary(ich_test_mod_sp)
# this is model 6 in the 'global' mod table btw
ich_from_tab_sp <-get.models(ichneumonidae_global_sp, subset = 6)[[1]]
summary(ich_from_tab_sp)


#### Geocoris ####
# Build global model selection table
geocoris_mod_table_sp <- readRDS('modTabs/Geocoris_spring_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = geocoris_mod_table_sp, rank = 1, 'Spring')
# Top 8 models are null. 9th doesn't look very good.
# Plot var importance charts
plotVarImportance(mod_table = geocoris_mod_table_sp)
# make global mod table (single distweight)
# try with null model
geocoris_global_sp <- makeGlobalLandcoverModTab(geocoris_mod_table_sp, 1)
# View(geocoris_global_sp)# top mod is null
buildBestLandcoverMod(geocoris_global_sp, 1,'Spring')
geocoris_shan_mod_sp <- buildFinalLandcoverMod(geocoris_mod_table_sp, geocoris_global_sp, 1)
summary(geocoris_shan_mod_sp)
# try with landcover vars
geocoris_global_sp <- makeGlobalLandcoverModTab(geocoris_mod_table_sp, 9)
buildBestLandcoverMod(geocoris_global_sp, 1)

png('spring_geo_cover_effect.jpg',
    width = 7,
    height = 5,
    units = 'in',
    res = 300)
geocoris_fin_mod_sp <- buildFinalLandcoverMod(geocoris_mod_table_sp, geocoris_global_sp, 6)
dev.off()
geocoris_fin_mod_sp <- buildFinalLandcoverMod(geocoris_mod_table_sp, geocoris_global_sp, 1)


summary(geocoris_fin_mod_sp)
plot(allEffects(geocoris_fin_mod_sp, residuals = TRUE), main = 'Geocoris, final spring model')
# best mod seems to be shan only

### fall data only ####

fD_fall <- field_data %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')


# define list of taxa of interest
taxaList <- c('Arachnida',
              'Coccinellidae',
              'Geocoris',
              'Ichneumonidae',
              'Acyrthosiphon',
              'NonAcy')

# read in tabs
tabList <- c()
for (i in 1:length(taxaList)) {

  tabList[[i]] <- readRDS(paste0('modTabs/', taxaList[[i]], '_fall_8class'))

}
# generate vi heatmaps
viList <- c()
for (i in 1:length(taxaList)) {

  viList[[i]] <- plotVarImportance(tabList[[i]], 'Fall')

}

viList[[1]]
viList[[2]]
viList[[3]]
viList[[4]]
viList[[5]]
viList[[6]]
# view head of tabs
heads <- c()
for (i in 1:length(taxaList)) {

  heads[[i]] <- head(tabList[[i]])

}

# plot top mod effects (one at a time, manually)

k =4
View(tabList[[1]])
plotName <- paste0(taxaList[[k]], '_rank1_mod_effect_fa.png')
plotMod <- get.models(tabList[[k]], 1)[[1]]
# plotMod <- get.models(acyG1, 4)[[1]]
summary(plotMod)
nullMod <- lmer(log(Ichneumonidae
                    + 1)~(1|Site), data = fD_fall, REML = F)
anova(plotMod, nullMod)
png(paste0('b',plotName),
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('ag_sig1',plotMod, residuals = TRUE),
     main = NULL,
     xlab = 'Agricultural cover,\n Very aggressive weighting',
     # ylim = c(0,10),
     id = list(n = length(get(plotMod@call$data)$id),
               labels = get(plotMod@call$data)$id))
dev.off()




ggplot(fD_fall, aes(reorder(id, Anthocoridae), log(Anthocoridae+1))) +
  geom_col()
  # theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=id), hjust = 0, vjust=0)+
  labs(x = 'Very aggressive', y = 'Moderate',
       title = 'Natural landcover, correlation between weightings')+
  coord_cartesian(xlim = c(-1,2.7))
ggsave('non-field.png', width = 6, height = 6, units = 'in')
# make global modsel tabs
araFD1 <- makeGlobalLandcoverModTabFall(tabList[[1]], 1, 3)
ggsave('anthG.png', width = 6, height = 6, units = 'in')
plotGlobalVarImportance(araFD1)
fD_fall_sub


#### Anthocoridae ####
buildLandcoverModTab <- function(taxon = 'empty', data = 'empty', m.max = 3){
  # taxon='NonAcy'
  # data=fD_fall

  distList <- c('_no ',
                '_const ',
                '_sig1 ',
                '_sig2 ',
                '_sig3 ',
                '_sig4 ',
                '_sig5 ')

  cand_mod_tabs <- list()

  if (taxon == 'empty' | !is_tibble(data)) {
    stop(red("Please specify taxon and data \n"), call. = FALSE)
  }

  # print inputs
  cat(yellow('Taxon:'),
      green(taxon),
      yellow('Data:'),
      green(deparse(substitute(data))),
      '\n')

  # get dfname
  dfname <- as.name(deparse(substitute(data)))
  # Build models
  for (i in 1:length(distList)) {
    # i=2

    # incase full model fails to fit, try 'tricking' dredge
    message(blue('fitting', distList[[i]],'models'))
    # fit reduced model
    fmod.red <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')
    # define full model formula
    form <-formula(
      paste0('log(',
             taxon,
             ' + 1) ~ ',
             paste0(varList, distList[[i]], '+ ', collapse = ''),
             '(1|Site)'
      ))
    # Replace reduced model formula with full global model formula.
    attr(fmod.red@frame, "formula") <- form
    # Run dredge() with m.max parameter to avoid convergence failures.
    fmod.red@call$data <- dfname

    cand_mod_tabs[[i]] <-  # superassign?
      model.sel(lapply(
        dredge(fmod.red, m.lim = c(NA, m.max), trace = 2, evaluate = FALSE),
        eval))

    message(blue(nrow(cand_mod_tabs[[i]]), 'models fit\n'))

  }
  # Rbind the elements of the list together. This forces recalculation of AICc
  mod_table <- rbind(cand_mod_tabs[[1]],
                     cand_mod_tabs[[2]],
                     cand_mod_tabs[[3]],
                     cand_mod_tabs[[4]],
                     cand_mod_tabs[[5]],
                     cand_mod_tabs[[6]],
                     cand_mod_tabs[[7]])

  return(mod_table)
}

# anthMods <- buildLandcoverModTab('Anthocoridae', fD_fall, 3)
plotName <- paste0('Anth', '_rank1_mod_effect_fa.png')
View(anthMods)
plotVarImportance(anthMods)
plotMod <- get.models(anthMods, 10)[[1]]
# plotMod <- get.models(acyG1, 4)[[1]]
summary(plotMod)
nullMod <- lmer(log(Anthocoridae
                    + 1)~(1|Site), data = fD_fall, REML = F)
anova(plotMod, nullMod)
png(paste0('a',plotName),
    width = 12.5,
    height = 6,
    units = 'in',
    res = 300)
plot(effect('alfalfa_no',plotMod, residuals = TRUE),
     main = NULL,
     xlab = 'Alfalfa cover,\n No weighting',
     id = list(n = length(get(plotMod@call$data)$id),
               labels = get(plotMod@call$data)$id))
dev.off()
#### Arachnida ####
# Build global model selection table
arachnida_mod_table_fa <- buildLandcoverModTab('Arachnida', fD_fall)
# NOT FITTING at all
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = arachnida_mod_table_fa, rank = 1, 'Fall')

# Plot var importance charts
plotVarImportance(mod_table = arachnida_mod_table_fa)
# make global mod table (single distweight)
arachnida_global_fa <- makeGlobalLandcoverModTab(arachnida_mod_table_fa, 1)
# plot var importance for global mod table

plotGlobalVarImportance(arachnida_global_fa)

# build best model
arachnida_fin_mod_fa <- buildFinalLandcoverMod(arachnida_mod_table_fa, arachnida_global_fa, 1) # Error
summary(arachnida_fin_mod_fa)
plot(allEffects(arachnida_fin_mod_fa, residuals = TRUE), main = 'Arachnida, final spring model')


#### Coccinellidae ####
# Build global model selection table
coccinellidae_mod_table_fa <- readRDS('modTabs/Coccinellidae_fall_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_fa, rank = 1,'Fall')
# Plot var importance charts
plotVarImportance(mod_table = coccinellidae_mod_table_fa)
# clear maximum at naturalArid_sig2
# check mod table
# View(coccinellidae_mod_table_fa)
# check sig2 top model
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_fa, rank = 2,'Fall')
# reminder - fall coccinellidae densities are mostly zero:
pred_data %>%
  filter(Taxa == 'Coccinellidae') %>%
  ggplot(aes(x = LogDensity, y = Site, fill = Season)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Coccinellidae density, Log+1 Transformation') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 6)
ggsave('fall_ladybug_logDens.jpg', width = 7, height = 5, units = 'in')
# maybe could do some zero-inf binomial models, but for now:

# make global mod table (single distweight)
coccinellidae_global_fa <- makeGlobalLandcoverModTabFall(coccinellidae_mod_table_fa, 1, 2)

# plot var importance for global mod table
png('fall_cocc_global_varImp.png',
    width = 7,
    height = 5,
    units = 'in',
    res = 300)
plotGlobalVarImportance(coccinellidae_global_fa)
dev.off()

test <- get.models(coccinellidae_global_fa, 1)[[1]]
test2 <- get.models(coccinellidae_global_fa, 2)[[1]]
summary(test)
summary(test2)
anova(test, test2)

# build "best" model
coccinellidae_fin_mod_fa <- buildFinalLandcoverMod(coccinellidae_mod_table_fa, coccinellidae_global_fa, 1)
summary(coccinellidae_fin_mod_fa)
plot(allEffects(coccinellidae_fin_mod_fa, residuals = TRUE),
     main = 'Coccinellidae, final spring model',
     id = list(n = length(get(coccinellidae_fin_mod_fa@call$data)$id),
               labels = get(coccinellidae_fin_mod_fa@call$data)$id))


test <- lmer(log(Coccinellidae+1)~naturalArid_sig2+(1|Site), data = fD_fall_sub, REML = F)
summary(test)
# singular
test2 <- lmer(log(Coccinellidae+1)~wet_sig2+(1|Site), data = fD_fall_sub, REML = F)
summary(test2)

mods <- list(test, test2)

testumod<-model.sel(mods)

# also singular, lower p value
dev.off()
plot(allEffects(test, residuals = TRUE),
     main = 'Coccinellidae, margin(no yerington) data model, fall',
     id = list(n = length(get(test@call$data)$id),
               labels = get(test@call$data)$id))


fD_fall %>% ggplot(aes(x = reorder(id, Coccinellidae), y = log(Coccinellidae + 1))) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Field')
ggsave(filename = 'coccFieldDensity.jpg', width = 7, height = 5, units = 'in')

fD_fall %>% ggplot(aes(x = reorder(id, Coccinellidae), y = naturalArid_sig2)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Field')
ggsave(filename = 'na4.jpg', width = 7, height = 5, units = 'in')

plot2mod <- lmer(log(Coccinellidae + 1)~wet_sig2+(1|Site), data = fD_fall, REML = TRUE)


## IDEA ####
# check greeness of different sD landcover classes at spring and fall dates, using sentinel imagery
# would give rough idea of which areas are greenest at different times.
ndviTable <- read_csv('raw_data/ndviStats.csv', col_types = 'ccc') %>%
  separate(ndviMean, c('meanSpring', 'meanFall'), sep = ',') %>%
  separate(ndviMedian, c('medianSpring', 'medianFall'), sep = ',') %>%
  separate(ndviStdDev, c('sdSpring', 'sdFall'), sep = ',') %>%
  mutate(across(.fns = ~
                  str_extract(., '[:digit:]\\.[:digit:]+E?-?[:digit:]?'))) %>%
  mutate(across(.fns = parse_number), .keep = 'none') %>%
  cbind(., klass = c('rye',
                     'bare',
                     'dairy',
                     'onion',
                     'paved',
                     'weedy',
                     'alfalfa',
                     'dirtRoad',
                     'riparian',
                     'structures',
                     'naturalArid',
                     'water')) %>%
  mutate(klass = as_factor(klass))

ggplot(ndviTable, aes(reorder(klass, desc(meanSpring)), meanSpring)) +
  geom_col() +
  ylim(c(0, 0.25)) +
  labs(x = 'Landcover class',
       y = 'Sentinel-2 NDVI, spring')
ggsave('ndvi.jpg', height = 5, width = 7, units = 'in')

ggplot(ndviTable, aes(reorder(klass, desc(meanSpring)), meanFall)) +
  geom_col() +
  ylim(c(0, 0.25))
ggplot(ndviTable, aes(reorder(klass, desc(meanSpring)), sdSpring)) +
  geom_col()
ggplot(ndviTable, aes(reorder(klass, desc(meanSpring)), sdFall)) +
  geom_col()

# weedy class is the greenest class

#### Ichneumonidae ####
# Build global model selection table
ichneumonidae_mod_table_fa <- buildLandcoverModTab('Ichneumonidae', fD_fall)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = ichneumonidae_mod_table_fa, rank = 1,'Fall')
# Plot var importance charts
plotVarImportance(mod_table = ichneumonidae_mod_table_fa)
# make global mod table (single distweight)
ichneumonidae_global_fa <- makeGlobalLandcoverModTab(ichneumonidae_mod_table_fa, 1)
# plot var importance for global mod table
plotGlobalVarImportance(ichneumonidae_global_fa)
# build best model
ichneumonidae_fin_mod_fa <- buildFinalLandcoverMod(ichneumonidae_mod_table_fa,
                                            ichneumonidae_global_fa,
                                            2)
summary(ichneumonidae_fin_mod_fa)
plot(allEffects(ichneumonidae_fin_mod_fa, residuals = TRUE), main = 'Ichneumonidae, final spring model')
# hypothetical ich model with dirt_const and shan:
ich_test_mod_fa <- lmer(log(Ichneumonidae + 1) ~ dirt_const + shan + (1|Site),
                        data = fD_fall_sub,
                        REML = FALSE)
summary(ich_test_mod_fa)
# this is model 6 in the 'global' mod table btw
ich_from_tab_fa <-get.models(ichneumonidae_global_fa, subset = 6)[[1]]
summary(ich_from_tab_fa)


#### Geocoris ####
# Build global model selection table
geocoris_mod_table_fa <- readRDS('modTabs/Geocoris_fall_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = geocoris_mod_table_fa, rank = 1, 'Fall')
# Top 8 models are null. 9th doesn't look very good.
# Plot var importance charts
plotVarImportance(mod_table = geocoris_mod_table_fa)
# make global mod table (single distweight)
# try with null model
geocoris_global_fa <- makeGlobalLandcoverModTab(geocoris_mod_table_fa, 1)
# View(geocoris_global_fa)# top mod is null
buildBestLandcoverMod(geocoris_global_fa, 1,'Fall')
geocoris_fin_mod_fa <- buildFinalLandcoverMod(geocoris_global_fa, geocoris_global_fa, 1)
summary(geocoris_fin_mod_fa)
# try with landcover vars
geocoris_global_fa <- makeGlobalLandcoverModTab(geocoris_mod_table_fa, 1)
buildBestLandcoverMod(geocoris_global_fa, 1)



### Aphids ~ landcover ####
#### Spring ####
##### Acyrthosiphon ####

# Build global model selection table
acyrthosiphon_mod_table_sp <- readRDS('modTabs/Acyrthosiphon_spring_8class')
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = acyrthosiphon_mod_table_sp, rank = 1, 'Spring')
# similar effect with sig1
buildBestLandcoverMod(mod_table = acyrthosiphon_mod_table_sp, rank = 1, 'Spring')
# Plot var importance charts
plotVarImportance(mod_table = acyrthosiphon_mod_table_sp)
# make global mod table (single distweight)
acyrthosiphon_global_sp <- makeGlobalLandcoverModTabSpring(acyrthosiphon_mod_table_sp, 1)
# plot var importance for global mod table
plotGlobalVarImportance(acyrthosiphon_global_sp)
# build best model
acyrthosiphon_fin_mod_sp <- buildFinalLandcoverMod(acyrthosiphon_mod_table_sp, acyrthosiphon_global_sp, 2)
summary(acyrthosiphon_fin_mod_sp)
plot(allEffects(acyrthosiphon_fin_mod_sp, residuals = TRUE),
     main = 'acyrthosiphon, final spring model',
     id = list(n = nrow(acyrthosiphon_fin_mod_sp@frame),
               labels = acyrthosiphon_fin_mod_sp@frame %>%
                 pull(Site)))
dev.off()
# View(acyrthosiphon_fin_mod_sp)

## add in coccinellidae
modSig5 <- lmer(log(Acyrthosiphon + 1) ~
                  alfalfa_sig5+naturalArid_sig5+
                  dirt_sig5+
                  ag_sig5+
                  impermeable_sig5+
                  weedy_sig5+wet_sig5+water_sig5+
                  log(Coccinellidae+1)+(1|Site),
                data = fD_spring,
                REML = FALSE,
                na.action = 'na.fail')
modSig3 <- lmer(log(Acyrthosiphon + 1) ~
                  alfalfa_sig3+naturalArid_sig3+
                  dirt_sig3+
                  ag_sig3+
                  impermeable_sig3+
                  weedy_sig3+wet_sig3+water_sig3+
                  log(Coccinellidae+1)+(1|Site),
                data = fD_spring,
                REML = FALSE,
                na.action = 'na.fail')
acyD1 <- dredge(modSig5, m.lim = c(1,4))
acyD2 <- dredge(modSig3, m.lim = c(1,4))
acyS5best <- get.models(acyD1, 1)[[1]]
acyS3best <- get.models(acyD2, 1)[[1]]
# acyS5bestPred <- get.models(acyD1pred, 1)[[1]]
View(acyD1)
View(acyD2)
plotGlobalVarImportance(acyD1)
plotGlobalVarImportance(acyD2)
ggsave('acyGlobal2.png', width = 12.5, height = 6, units = 'in')
acyMod <- lmer(log(Acyrthosiphon + 1) ~

                 log(Coccinellidae+1)+(1|Site),
               data = fD_spring,
               REML = FALSE,
               na.action = 'na.fail')
summary(acyMod)

#### Fall ####
##### Non-acyrthosiphon #####
# Build global model selection table
nonacy_mod_table_fl <- buildLandcoverModTab('NonAcy', fD_fall)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = nonacy_mod_table_fl, rank = 1, 'Spring')
# similar effect with sig1
buildBestLandcoverMod(mod_table = nonacy_mod_table_fl, rank = 8, 'Spring')
# Plot var importance charts
plotVarImportance(mod_table = nonacy_mod_table_fl)
# make global mod table (single distweight)
nonacy_global_fl <- makeGlobalLandcoverModTab(nonacy_mod_table_fl, 1)
# plot var importance for global mod table
plotGlobalVarImportance(nonacy_global_fl)
# build best model
nonacy_fin_mod_fl <- buildFinalLandcoverMod(nonacy_mod_table_fl, nonacy_global_fl, 2)
summary(nonacy_fin_mod_fl)
plot(allEffects(nonacy_fin_mod_fl, residuals = TRUE),
     main = 'nonacy, final spring model',
     id = list(n = nrow(nonacy_fin_mod_fl@frame),
               labels = nonacy_fin_mod_fl@frame %>%
                 pull(Site)))
# View(nonacy_fin_mod_fl)

# Spring-Fall deltas ####

# create delta columns
delta_density <- mean_density_field %>%
  pivot_wider(names_from = Season, values_from = 5:18) %>%
  filter(Site != 'Minden') %>%
  mutate(Acyrthosiphon = Acyrthosiphon_Fall - Acyrthosiphon_Spring,
         Therioaphis = Therioaphis_Fall - Therioaphis_Spring,
         Aphis = Aphis_Fall - Aphis_Spring,
         AllAph = AllAph_Fall - AllAph_Spring,
         NonAcy = NonAcy_Fall - NonAcy_Spring,
         Anthocoridae = Anthocoridae_Fall - Anthocoridae_Spring,
         Arachnida = Arachnida_Fall - Arachnida_Spring,
         Coccinellidae = Coccinellidae_Fall - Coccinellidae_Spring,
         Geocoris = Geocoris_Fall - Geocoris_Spring,
         Ichneumonidae = Ichneumonidae_Fall - Ichneumonidae_Spring,
         Nabis = Nabis_Fall - Nabis_Spring) %>%
  select(id, Site, Field, Acyrthosiphon, Therioaphis, Aphis, AllAph, NonAcy,
         Anthocoridae, Arachnida, Coccinellidae, Geocoris, Ichneumonidae,
         Nabis)

delta.mod <- lmer(AllAph ~ Anthocoridae + Arachnida + Coccinellidae +
                    Geocoris + Ichneumonidae + Nabis + (1|Site),
                  data = delta_density,
                  na.action = 'na.fail',
                  REML = FALSE)
summary(delta.mod)
delta_dredge <- dredge(delta.mod, m.max = 3)

# Plot the variable importance.
importance_tab <- sw(delta_dredge) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    arrange(sw)

p <- ggplot(data = importance_tab, aes(x = reorder(names, -sw), y = '', fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Predator',
       y = '',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Extract and examine the best model.
# Best mod is null. 2nd best below.
best.delta.mod <- get.models(delta_dredge, subset = 2)[[1]]
plot(allEffects(best.delta.mod, residuals = TRUE))

# no good!! note that AICc delta is 6.69 from best (null) mod.
summary(best.delta.mod)

# SEM ####
## Sopring ####
### Acyrthosiphon ####
# review models
summary(best.acy.mod.sp)
summary(acyrthosiphon_fin_mod_sp)
# need coccinellidae, richness?
summary(coccinellidae_fin_mod_sp)
# need weedyWet_sig4
# format data
# need: plot-level insect densities, no margin data
landcover_sem <- landcover_wide %>%
  mutate(field = substr(field, 2, 2),
         fieldID = paste(site, field))

sem_data <- data_long %>%
  mutate(Density = log(Density + 1)) %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  mutate(fieldID = paste(Site, Field)) %>%
  left_join(landcover_sem, by = 'fieldID') %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-') %>%
  select(-`id.x`) %>%
  mutate(measurementID = paste(Site, Field, Plot)) %>%
  pivot_wider(              names_from = c(Taxa),
              values_from = Density) %>%
  cbind(., to.dummy(.$Treatment, 'Trt'))

sem_data_fall <- data_long %>%
  mutate(Density = log(Density + 1)) %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  mutate(fieldID = paste(Site, Field)) %>%
  left_join(landcover_sem, by = 'fieldID') %>%
  filter(Season == 'Fall',
         Treatment != 'Pre-') %>%
  select(-`id.x`) %>%
  mutate(measurementID = paste(Site, Field, Plot)) %>%
  pivot_wider(              names_from = c(Taxa),
                            values_from = Density) %>%
  cbind(., to.dummy(.$Treatment, 'Trt'))


names(sem_data)
fD_spring_flood$id

sem_data2 <- sem_data %>% left_join(fD_spring_flood %>% select(id, wateringMethod), by = c('id.y'='id')) %>%
  cbind(., to.dummy(.$wateringMethod, 'Wtr'))


detach("package:lmerTest", unload=TRUE) # this fucks with psem for some reason
spring_acy_sem <- psem(
  lmer(Coccinellidae ~ `weedy_sig1` + `Trt.Sham` + Acyrthosiphon +(1|site),
     data = sem_data2),
  lmer(Acyrthosiphon ~ `Wtr.Flooding` +(1|site),
     data = sem_data2)
  )
plot(spring_acy_sem)
qinstall.packages('DiagrammeRsvg')

plot(spring_acy_sem)%>%
  export_svg %>% charToRaw %>% rsvg_png("graph.png")
install.packages('rsvg')

export_graph(plot(spring_acy_sem),
             file_name = "pic.png",
             file_type = "png")
plot(spring_acy_sem, show = 'unstd') #what do standardized coefs mean?


fall_acy_sem <- psem(
  lmer(Acyrthosiphon ~ Ichneumonidae + `Trt.Sham` + (1|Site),
       data = sem_data_fall),
  # lmer(Anthocoridae ~ `ag_sig1` + `naturalArid_sig1` + `Trt.Sham` + (1|Site),
  #      data = sem_data_fall),
  lmer(Ichneumonidae ~ `dirt_sig4` + `water_sig4` + `Trt.Sham` + (1|Site),
       data = sem_data_fall)
)


fall_nonacy_sem <- psem(
  lmer(NonAcy ~ Arachnida + Ichneumonidae + `Trt.Sham` + (1|Site),
       data = sem_data_fall),
  lmer(Arachnida ~ `ag_sig2` + `Trt.Sham` + (1|Site),
       data = sem_data_fall),
  lmer(Ichneumonidae ~ `dirt_sig4` + `water_sig4` + `Trt.Sham` + (1|Site),
       data = sem_data_fall)
)

foo <-plot(fall_acy_sem)
bar <-plot(fall_nonacy_sem)

foo
bar
foo %>%
export_svg() %>%
  charToRaw %>%
  rsvg_png("./foo.svg")
bar %>%
export_svg() %>%
  charToRaw %>%
  rsvg_png("./bar.svg")

field_data %>%
  select(id, Season, Coccinellidae) %>%
  pivot_wider(id_cols = id, names_from = Season,
              values_from = Coccinellidae) %>%
  ggplot(aes(Spring, Fall, label = id)) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0) +
  geom_smooth(method = 'lm')#+
  # geom_text(label = 'R^2 = 0.09')
  # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)
ggsave('s_f_c.png', width = 12.5, height = 6, units = 'in')



plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, "7class-fixed-savedplots/", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("7class-fixed-savedplots/", 1:length(sorted.png.names), ".png")

# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
file.copy(from=plots.png.paths, to="7class-fixed-savedplots")
file.rename(from=sorted.png.names, to=numbered.png.names)


# functions ####
buildLandcoverModTab <- function(taxon = 'empty', data = 'empty', m.max = 5){
  # taxon='NonAcy'
  # data=fD_fall

  distList <- c('_no ',
                '_const ',
                '_sig1 ',
                '_sig2 ',
                '_sig3 ',
                '_sig4 ',
                '_sig5 ')

  cand_mod_tabs <- list()

  if (taxon == 'empty' | !is_tibble(data)) {
    stop(red("Please specify taxon and data \n"), call. = FALSE)
  }

  # print inputs
  cat(yellow('Taxon:'),
      green(taxon),
      yellow('Data:'),
      green(deparse(substitute(data))),
      '\n')

  # get dfname
  dfname <- as.name(deparse(substitute(data)))
  # Build models
  for (i in 1:length(distList)) {
    # i=2

    # incase full model fails to fit, try 'tricking' dredge
    message(blue('fitting', distList[[i]],'models'))
    # fit reduced model
    fmod.red <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')
    # define full model formula
    form <-formula(
      paste0('log(',
             taxon,
             ' + 1) ~ ',
             paste0(varList, distList[[i]], '+ ', collapse = ''),
             '(1|Site)'
      ))
    # Replace reduced model formula with full global model formula.
    attr(fmod.red@frame, "formula") <- form
    # Run dredge() with m.max parameter to avoid convergence failures.
    fmod.red@call$data <- dfname

    cand_mod_tabs[[i]] <-  # superassign?
      model.sel(lapply(
        dredge(fmod.red, m.lim = c(NA, m.max), trace = 2, evaluate = FALSE),
        eval))

    message(blue(nrow(cand_mod_tabs[[i]]), 'models fit\n'))

  }
  # Rbind the elements of the list together. This forces recalculation of AICc
  mod_table <- rbind(cand_mod_tabs[[1]],
                     cand_mod_tabs[[2]],
                     cand_mod_tabs[[3]],
                     cand_mod_tabs[[4]],
                     cand_mod_tabs[[5]],
                     cand_mod_tabs[[6]],
                     cand_mod_tabs[[7]])

  return(mod_table)
}

#### Anthocoridae spring ####
# Build global model selection table
Anthocoridae_mod_table_sp <- buildLandcoverModTab('Anthocoridae', fD_spring)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = Anthocoridae_mod_table_sp, rank = 1, 'Spring')

# Plot var importance charts
plotVarImportance(mod_table = Anthocoridae_mod_table_sp)
# make global mod table (single distweight)
Anthocoridae_global_sp <- makeGlobalLandcoverModTab(Anthocoridae_mod_table_sp, 1)
# View(Anthocoridae_global_sp)# top mod is null
buildBestLandcoverMod(Anthocoridae_global_sp, 1,'Spring')
Anthocoridae_fin_mod_sp <- buildFinalLandcoverMod(Anthocoridae_global_sp, Anthocoridae_global_sp, 1)
summary(Anthocoridae_fin_mod_sp)

test <- lmer(log(Anthocoridae + 1)~total_cover+(1|Site), data = fD_spring, REML = FALSE)

summary(test)
plot(allEffects(test, residuals = TRUE))


stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,

                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))

                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }

                            params
                          },

                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }

                            if (is.null(data$weight)) data$weight <- 1

                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }

                            if (is.character(method)) method <- match.fun(method)

                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))

                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                             list(a = format(coef(m)[1], digits = 3),
                                                  b = format(coef(m)[2], digits = 3),
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))

                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)

                          },

                          required_aes = c("x", "y")
)




