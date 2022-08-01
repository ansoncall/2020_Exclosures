# Load packages ####
library(car) # for Anova() on lmer model objects
library(crayon) # for colored terminal outputs
library(effects) # for effects plots
library(emmeans) # for computing SEM marginal means
library(gridExtra) # create multi-panel plots
library(ggridges) # for ggridges plots
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
library(sjPlot) # create effects plots on lmer objects
library(tidyverse) # R packages for data science
library(varhandle) # easily create dummy vars with to.dummy()
library(vegan) # for diversity indices in vegdata

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

# read in data ####
data <- read_csv('tidy_data/data.csv', col_types = 'fffffdddddddddddddddddddd')
data_long <- read_csv('tidy_data/data_long.csv', col_types = 'fffffffd')
landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffdddddddddddd')
vegPlots <- read_csv('tidy_data/vegPlots.csv', col_types = 'fffffff')
vegSites <- read_csv('tidy_data/vegSites.csv', col_types = 'f')


# make df of total abundance across treatments
mean_density <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  separate(id, c('Site', 'Field', 'Plot', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density)

# define convenience identifiers ####
# define lists of predator and aphid taxa
predlist <- c('Arachnida','Coccinellidae','Ichneumonidae',
              'Nabis', 'Geocoris', 'Anthocoridae')
aphlist <- c('Acyrthosiphon', 'Aphis', 'Therioaphis', 'AllAph', 'NonAcy')

# Arthropod data summary ####
## Aphid histograms ####
aph_data <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  filter(Taxa %in% aphlist) %>%
  mutate(LogDensity = log(Density+1))

# build density plot
ggplot(data = aph_data, aes(y = Taxa, x = LogDensity, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Aphids, Log+1 Transformation',
       x = 'log(Density)',
       y = 'Taxon') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 10)

# barchart for slideshow
# spring
ggplot(data = aph_data %>%
         filter(Season == 'Spring',
                Taxa %in% c('Acyrthosiphon', 'Aphis', 'Therioaphis')),
       aes(x = LogDensity, y = Taxa, fill = Taxa)) +
  geom_boxplot() +
  xlim(c(0, 8.5)) +
  scale_fill_manual(values = c('#548235', '#3b3838', '#c78f00'),
                    guide = 'none') +
  labs(y = 'Aphid genus', x = 'log(Density)')+
  coord_flip() +
  theme(text = element_text(size = 20))
# ggsave('spring_aphid_density.jpg', width = 7, height = 5)
# fall
ggplot(data = aph_data %>%
         filter(Season == 'Fall',
                Taxa %in% c('Acyrthosiphon', 'Aphis', 'Therioaphis')),
       aes(x = LogDensity, y = Taxa, fill = Taxa)) +
  geom_boxplot() +
  xlim(c(0, 8.5)) +
  scale_fill_manual(values = c('#548235', '#3b3838', '#c78f00'),
                    guide = 'none') +
  labs(y = 'Aphid genus', x = 'log(Density)') +
  coord_flip() +
  theme(text = element_text(size = 20))
# ggsave('fall_aphid_density.jpg', width = 7, height = 5)

# total
ggplot(data = aph_data %>%
         filter(Taxa == 'AllAph'),
       aes(x = LogDensity, y = Season, fill = Season)) +
  geom_boxplot() +
  xlim(c(0, 8.5)) +
  scale_fill_manual(values = c('#008000','#993300'),
                    guide = 'none') +
  labs(y = 'Aphid genus', x = 'log(Density)') +
  coord_flip() +
  theme(text = element_text(size = 20),
        axis.text.x=element_text(angle=30,hjust=0.9)) +
  labs(x = element_blank(), y = element_blank())
# ggsave('total_aphid_density.jpg', width = 1.5, height = 4.5)

## Predator histogram ####
pred_data <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  filter(Taxa %in% predlist) %>%
  mutate(LogDensity = log(Density+1), Taxa = fct_relevel(Taxa, sort))

# build density plot
ggplot(data = pred_data, aes(y = Taxa, x = LogDensity, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Predators, Log+1 Transformation') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 6)

# barchart (density*season) for slideshow
ggplot(data = pred_data,
       aes(x = LogDensity, y = Taxa, fill = Season)) +
  geom_boxplot() +
  coord_flip() +
  xlim(c(0, 6)) +
  labs(y = 'Predator taxon', x = 'log(Density)') +
  scale_fill_manual(values = c('#008000','#993300')) +
  theme(text = element_text(size = 20),
        axis.text.x=element_text(angle=30,hjust=0.9))
# ggsave('spring_pred_density.jpg', width = 10, height = 5)

# demo figure for field effect on density
ggplot(data = pred_data %>%
         filter(Taxa == 'Coccinellidae',
                Site == 'Minden',
                Field %in% c(1, 2)),
       aes(x = LogDensity, y = Field, fill = Taxa)) +
  geom_boxplot() +
  coord_flip() +
  guides(fill = 'none') +
  labs(y = 'Field', x = 'log(Density)', title = 'Ladybug density') +
  theme(text = element_text(size = 20))
# ggsave('example_field_effect.jpg', width = 4, height = 5)

# Landcover data summary ####

# define varlist
varList <- c('alfalfa',
             'naturalArid',
             'dirt',
             'ag',
             'impermeable',
             'weedy',
             'wet',
             'water')

# Collapse classes
landcover %<>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedy = class5,
         wet = class8,
         water = class11
         ) %>%
  select(-(class0:class11)) %>%
  mutate(distanceWeight = fct_relevel(distanceWeight, 'no', after = Inf))
levels(landcover$distanceWeight)
# lengthen
landcover_long <- landcover %>%
  pivot_longer(alfalfa:water,
               names_to = 'class',
               values_to = 'areaScore') %>%
  select(siteId, distanceWeight, site, field, class, areaScore)

ggplot(landcover_long,
       # aes - limit x axis to 4 characters
       aes(x = sub(class, pattern = "(\\w{4}).*", replacement = "\\1."),
           y = areaScore,
           fill = class)) +
  geom_bar( stat = "identity") +
  facet_wrap(ncol = 7, ~site*distanceWeight, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Notes #### Obviously this can't work. Need a way to succinctly show how
## distance weighting affects the areaScores. Also, color, class names, etc.
## need to be fixed. Messaged Beth for advice.

# Margin data summary ####
## Across sites ####
# Shannon diversity
ggplot(data = vegPlots %>% filter(type == 'Margin'),
       aes(x = site, y = shan, fill = season)) +
  geom_boxplot() +
  labs(title = 'Shannon diversity index',
       x = 'Site', y = 'Diversity')
# Plant species richness
ggplot(data = vegPlots %>% filter(type == 'Margin'),
       aes(x = site, y = rich, fill = season)) +
  geom_boxplot() +
  labs(title = 'Plant species richness',
       x = 'Site', y = 'Richness')
# Total plant cover
ggplot(data = vegPlots %>% filter(type == 'Margin') %>%
         mutate(total_cover = select(., 12:132) %>% rowSums(na.rm = TRUE)),
       aes(x = site, y = total_cover, fill = season)) +
  geom_boxplot() +
  labs(title = 'Plant cover %',
       x = 'Site', y = '%')
### Notes ####
## Not seeing much variation here. May want to shaw by field.
## Remember, Yerington will always be missing. Need to fix colors, factor order,
## etc.

# Model selection ####
## Aphids ~ predators ####
### spring data only ####

# define list of aphid taxa
short_aphlist <- c('Acyrthosiphon',
             'Aphis',
             'Therioaphis')

# create empty lists
dredges <- list()

for (i in 1:length(short_aphlist)) {

  taxon <- short_aphlist[[i]]
  formula <-  paste0('log(',
                     taxon,
                     ' + 1) ~ log(Arachnida + 1) + log(Coccinellidae + 1) + ',
                     'log(Ichneumonidae + 1) + log(Nabis + 1) + log(Geocoris + 1) +',
                     'log(Anthocoridae + 1) + (1|Site) + (1|Field)')
  mGlobal <- mean_density %>%
    filter(Season=='Spring') %$%
    lmer(formula(formula),
         na.action = 'na.fail', # can't replicate previous error?
         REML = FALSE)

  dredges[[i]] <- dredge(mGlobal) # dredge and store output in list

}
names(dredges) <- short_aphlist

# check data sources
importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(ExpVars = fct_recode(names,
                                Anthocoridae = 'log(Anthocoridae + 1)',
                                Arachnida = 'log(Arachnida + 1)',
                                Coccinellidae = 'log(Coccinellidae + 1)',
                                Geocoris = 'log(Geocoris + 1)',
                                Ichneumonoidea = 'log(Ichneumonidae + 1)',
                                Nabis = 'log(Nabis + 1)'),
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

# extract and examine the top Acyrthosiphon model.
best.acy.mod <- get.models(springTabs[[1]], subset = 1)[[1]]
summary(best.acy.mod)
# must remake to plot effects.
best.acy.mod.sp <- lmer(log(Acyrthosiphon + 1) ~ log(Coccinellidae + 1) +
                                                     (1|Site) +
                                                     (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'),
                     REML = FALSE, na.action = 'na.omit')
summary(best.acy.mod.sp)
# Note: no data seems to be missing here.

# make plot
# png('spring_acy_effect.jpg',
#     width = 7,
#     height = 5,
#     units = 'in',
#     res = 300)
plot(allEffects(best.acy.mod.sp, residuals = TRUE),
     main = 'Acyrthosiphon, Spring',
     id = list(n = 36, labels = mean_density %>%
                 filter(Season == 'Spring') %>%
                 unite('id', Site, Field, Plot, sep = '.') %>%
                 pull(id)))
# dev.off()

# extract and examine the top Aphis model.
best.aphis.mod <- get.models(springTabs[[2]], subset = 1)[[1]]
summary(best.aphis.mod)
# Must remake to plot effects.
best.aphis.mod <- lmer(log(Aphis + 1) ~ log(Ichneumonidae + 1) +
                         (1|Site) +
                         (1|Field),
                       data = mean_density %>% filter(Season == 'Spring'),
                       REML = FALSE)
# make plot
# png('spring_aphis_effect.jpg',
#     width = 7,
#     height = 5,
#     units = 'in',
#     res = 300)
plot(allEffects(best.aphis.mod, residuals = TRUE),
     main = 'Aphis, Spring',
     id = list(n = 36, labels = mean_density %>%
                 filter(Season == 'Spring') %>%
                 unite('id', Site, Field, Plot, sep = '.') %>%
                 pull(id)))
# dev.off()

# extract and examine the top Therioaphis model.
best.therio.mod <- get.models(springTabs[[3]], subset = 1)[[1]]
summary(best.therio.mod)
# Must remake to plot effects.
best.therio.mod <- lmer(log(Therioaphis + 1) ~ log(Geocoris + 1) +
                          (1|Site) +
                          (1|Field),
                        data = mean_density %>% filter(Season == 'Spring'),
                        REML = FALSE)
# make plot
# png('spring_therio_effect.jpg',
#     width = 7,
#     height = 5,
#     units = 'in',
#     res = 300)
plot(allEffects(best.therio.mod, residuals = TRUE),
     main = 'Therioaphis, Spring',
     id = list(n = 36, labels = mean_density %>%
                 filter(Season == 'Spring') %>%
                 unite('id', Site, Field, Plot, sep = '.') %>%
                 pull(id)))
# dev.off()

# make variable importance heatmap
p <- bind_rows(importance_tab, .id = 'Taxon') %>%
  mutate(VarWeights = as.numeric(VarWeights),
         ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
         Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
  filter(ExpVars %in% c("Coccinellidae", "Geocoris", "Ichneumonoidea")) %>%
  ggplot(aes(x = ExpVars, y = Taxon, fill = VarWeights)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        legend.title = element_text('Variable Importance')) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = '\"Beneficial\" taxon', y = 'Aphid genus', title = 'Spring') +
  theme(text = element_text(size = 15),
        axis.text.x=element_text(angle=30,hjust=0.9))
p
# ggsave('spring_aphid_varweights.jpg', width = 6, height = 4)

### fall data only ####

# create empty lists
dredges <- list()

for (i in 1:length(short_aphlist)) {

  aph_density <- as.name(short_aphlist[[i]])
  mGlobal <- mean_density %>%
    filter(Season=='Fall') %$%
    lmer(log(eval(aph_density) + 1) ~ (log(Arachnida + 1) +
                                         log(Coccinellidae + 1) +
                                         log(Ichneumonidae + 1) +
                                         log(Nabis + 1) +
                                         log(Geocoris + 1) +
                                         log(Anthocoridae + 1) +
                                         (1|Site) +
                                         (1|Field)),
         na.action = 'na.fail',
         REML = FALSE)
  dredges[[i]] <- dredge(mGlobal) # dredge and store output in list

}
names(dredges) <- short_aphlist

importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tibble(names = names(.),
           .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(ExpVars = fct_recode(names,
                                Anthocoridae = 'log(Anthocoridae + 1)',
                                Arachnida = 'log(Arachnida + 1)',
                                Coccinellidae = 'log(Coccinellidae + 1)',
                                Geocoris = 'log(Geocoris + 1)',
                                Ichneumonoidea = 'log(Ichneumonidae + 1)',
                                Nabis = 'log(Nabis + 1)'),
           VarWeights = sw,
           .keep = 'none') %>%
    arrange(ExpVars)
})
names(importance_tab) <- short_aphlist
# make tables of model selection results
fallTabs <- lapply(dredges, slice_head, n = 5)
names(fallTabs) <- short_aphlist
# show tables
# fallTabs

# Extract and examine the top Acyrthosiphon model.
best.acy.mod <- get.models(fallTabs[[1]], subset = 1)[[1]]
summary(best.acy.mod)
# Must remake to plot effects.
best.acy.mod <- lmer(log(Acyrthosiphon + 1) ~ log(Anthocoridae + 1) +
                       log(Ichneumonidae + 1) +
                       (1|Site) +
                       (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'),
                     REML = FALSE)
# make plot
# png('fall_acy_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Fall')
# dev.off()

# Extract and examine the top Aphis model.
best.aphis.mod <- get.models(fallTabs[[2]], subset = 1)[[1]]
summary(best.aphis.mod)
# Must remake to plot effects.
best.aphis.mod <- lmer(log(Aphis + 1) ~ log(Arachnida + 1) +
                       log(Ichneumonidae + 1) +
                       (1|Site) +
                       (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'),
                     REML = FALSE)
# make plot
# png('fall_aphis_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
plot(allEffects(best.aphis.mod, residuals = TRUE), main = 'Aphis, Fall')
# dev.off()

# Extract and examine the top Therioaphis model.
best.therio.mod <- get.models(fallTabs[[3]], subset = 1)[[1]]
summary(best.therio.mod)
# Must remake to plot effects.
best.therio.mod <- lmer(log(Therioaphis + 1) ~ log(Arachnida + 1) +
                       (1|Site) +
                       (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'),
                     REML = FALSE)
# make plot
# png('fall_therio_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
plot(allEffects(best.therio.mod,
                residuals = TRUE),
     main = 'Therioaphis, Fall')
# dev.off()

# make plot
p <- bind_rows(importance_tab, .id = 'Taxon') %>%
  mutate(VarWeights = as.numeric(VarWeights),
         ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
         Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
  filter(ExpVars %in% c("Arachnida", "Geocoris", "Ichneumonoidea")) %>%
  ggplot(aes(x = ExpVars, y = Taxon, fill = VarWeights)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        legend.title = element_text('Variable Importance')) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = '\"Beneficial\" taxon', y = 'Aphid genus', title = 'Fall') +
  theme(text = element_text(size = 15),
        axis.text.x=element_text(angle=30,hjust=0.9))

p
# ggsave('fall_aphid_varweights.jpg', width = 6, height = 4)

# clean environment
rm(dredges, fallTabs,
   importance_tab, mGlobal, springTabs, aph_density)

# Predators ~ ####
### Prepare data ####
# summarize mean_density across fields
mean_density_field <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  separate(id, c('Site', 'Field', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density) %>%
  mutate(id = paste0(Site, 0, Field), .before = Site)
# Reshape landcover data
landcover_wide <- landcover %>%
  pivot_wider(id_cols = c(site, field),
              names_from = distanceWeight,
              values_from = alfalfa:water) %>%
  mutate(id = paste0(site, field), .before = site)
# join landcover data
field_data <- left_join(mean_density_field, landcover_wide)

# create margin data
field_margins <- vegPlots %>% filter(type == 'Margin') %>%
  mutate(total_cover = select(., 12:132) %>% rowSums(na.rm = TRUE)) %>%
  group_by(field_id, season) %>%
  summarize(shan = mean(shan),
            rich = mean(rich),
            total_cover = mean(total_cover),
            .groups = 'keep') %>%
  # make join key
  mutate(id = str_replace(field_id, ' ', '0'), Season = season) %>%
  ungroup() %>%
  select(-field_id, -season)

# join margin data
field_data <- left_join(field_data, field_margins)

# clean environment
rm(field_margins)

#### formulas ####
# formula for calculating model selection table for predators
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
  plot(allEffects(best.mod, residuals = TRUE),
       main = paste(paste(getPred, collapse = ''),
                    'landcover model rank',
                    as.character(rank),
                    'Season:', season),
       id = list(n = length(get(best.mod@call$data)$id),
                 labels = get(best.mod@call$data)$id))
  summary(best.mod)
}

# plot the variable importance
plotVarImportance <- function (mod_table, season = 'unknown season'){

  getPred <- get.models(mod_table, subset = 1)[[1]]@call$formula[[2]]
  importance_tab <- sw(mod_table) %>% #tibble(names = names(.))
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    arrange(names) %>%
    separate(names, c('class', 'distWeight'), sep = "_") %>%
    mutate(distWeight = as_factor(recode(distWeight,
                                         `no` = 'no decay',
                                         `const` = 'constant',
                                         `sig1` = 'very aggressive',
                                         `sig2` = 'aggressive',
                                         `sig3` = 'moderately aggressive',
                                         `sig4` = 'moderate',
                                         `sig5` = 'slight',
                                         `sig6` = 'minimal')))

  p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
    geom_tile() +
    theme(axis.text.x=element_text(angle = 45, hjust = 0),
          axis.text.y = element_text(angle = 45)) +
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         title = paste(paste(getPred, collapse = ''),
                       'variable importance',
                       season))

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
         fill = 'Variable importance')

  qq <- ggplotly(q, tooltip = 'weight')
  subplot(pp, qq, widths = c(7/8, 1/8))

}

# Add models with margin data to model selection table.
makeGlobalLandcoverModTab <- function(mod_table, rank = 1){
  mod_table = arachnida_mod_table_sp
  rank = 1
  mod <- get.models(mod_table, subset = rank)[[1]]
  # Global model is rank-deficient. We will have to 'trick' dredge.

  # mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
  #                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
  #                     (1|Site),
  #                  data = fD_spring_sub) # This doesn't work.
  response <- names(mod@frame) %>% head(1)
  landCoverPredictors <- names(mod@frame) %>% head(-1) %>% tail(-1)
  distWeight <- str_match(landCoverPredictors[[1]], '_[:alnum:]+')
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
  form.red <- formula('log(Arachnida + 1) ~ (1|Site)')
  # Fit reduced model.
  fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE,
                   na.action = 'na.fail')
  # Replace reduced model formula with full global model formula.
  attr(fmod.red@frame, "formula") <- form
  # Check formula attribute of reduced model.
  formula(fmod.red) # Looks good.

  # Run dredge() with m.max parameter to avoid convergence failures.
  ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
  ms1
}

# Plot the variable importance.
plotGlobalVarImportance <- function(global_tab){
  importance_tab <- sw(global_tab) %>%
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    mutate(names = str_replace(names, '[:digit:]', '')) %>%
    rename(varWeight = sw)

  p <- ggplot(data = importance_tab,
              aes(x = reorder(names,-varWeight), y = varWeight)) +
    geom_col() +
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
    labs(x = 'Variable', y = 'Importance')

  ggplotly(p, tooltip = 'sw')
}

# Build final model based on variable selection exercise.
buildFinalLandcoverMod <- function(no_margin_tab, global_tab, rank = 1){

  gMod <- get.models(global_tab, subset = rank)[[1]]
  globalPredictors <- names(gMod@frame) %>% head(-1) %>% tail(-1)
  nmMod <- get.models(no_margin_tab, subset = rank)[[1]]
  noMarginPredictors <- names(nmMod@frame) %>% head(-1) %>% tail(-1)
  response <- names(gMod@frame) %>% head(1)
  if (length(noMarginPredictors) == 0 && length(globalPredictors) == 0){
    cat(red('Top models are null, with or without margin data!'))
  } else if (length(globalPredictors) == 0){
    cat(red('Warning: Best global model is null\nFitting without field margin data for extra power \n'))
    cat(yellow('Using full dataset'))
    fin.mod <- lmer(as.formula(paste0(response, '~', noMarginPredictors,'+ (1|Site)')),
                    data = fD_spring, REML = FALSE)
  } else {
    cat(yellow('Margin data included in model;
               Dropping Yerington data'))
    fin.mod <- lmer(as.formula(paste0(response, '~', globalPredictors,'+ (1|Site)')),
                    data = fD_spring_sub, REML = FALSE)
  }
  fin.mod

}

### spring data only ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')

#### Arachnida ####
# Build global model selection table
arachnida_mod_table_sp <- buildLandcoverModTab('Arachnida', fD_spring)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = arachnida_mod_table_sp, rank = 1, 'Spring')
# Plot var importance charts
plotVarImportance(mod_table = arachnida_mod_table_sp)
# make global mod table (single distweight)
arachnida_global_sp <- makeGlobalLandcoverModTab(arachnida_mod_table_sp, 1)
# plot var importance for global mod table
plotGlobalVarImportance(arachnida_global_sp)
# build best model
arachnida_fin_mod_sp <- buildFinalLandcoverMod(arachnida_mod_table_sp, arachnida_global_sp, 2) # Error
summary(arachnida_fin_mod_sp)
plot(allEffects(arachnida_fin_mod_sp, residuals = TRUE), main = 'Arachnida, final spring model')


#### Coccinellidae ####
# Build global model selection table
coccinellidae_mod_table_sp <- buildLandcoverModTab('Coccinellidae', fD_spring)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_sp, rank = 1,'Spring')
# Plot var importance charts
plotVarImportance(mod_table = coccinellidae_mod_table_sp)
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
ichneumonidae_mod_table_sp <- buildLandcoverModTab('Ichneumonidae', fD_spring)
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
geocoris_mod_table_sp <- buildLandcoverModTab('Geocoris', fD_spring)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = geocoris_mod_table_sp, rank = 9, 'Spring')
# Top 8 models are null. 9th doesn't look very good.
# Plot var importance charts
plotVarImportance(mod_table = geocoris_mod_table_sp)
# make global mod table (single distweight)
# try with null model
geocoris_global_sp <- makeGlobalLandcoverModTab(geocoris_mod_table_sp, 1)
# View(geocoris_global_sp)# top mod is null
buildBestLandcoverMod(geocoris_global_sp, 2,'Spring')
geocoris_shan_mod_sp <- buildFinalLandcoverMod(geocoris_global_sp, geocoris_global_sp, 2)
summary(geocoris_shan_mod_sp)
# try with landcover vars
geocoris_global_sp <- makeGlobalLandcoverModTab(geocoris_mod_table_sp, 9)
buildBestLandcoverMod(geocoris_global_sp, 2)
geocoris_fin_mod_sp <- buildFinalLandcoverMod(geocoris_mod_table_sp, geocoris_global_sp, 2)
summary(geocoris_fin_mod_sp)
plot(allEffects(geocoris_fin_mod_sp, residuals = TRUE), main = 'Geocoris, final spring model')
# best mod seems to be shan only

### fall data only ####

fD_fall <- field_data %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')

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
coccinellidae_mod_table_fa <- buildLandcoverModTab('Coccinellidae', fD_fall)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_fa, rank = 1,'Fall')
# Plot var importance charts
plotVarImportance(mod_table = coccinellidae_mod_table_fa)
# extremely bipolar distribution??
# check mod table
# View(coccinellidae_mod_table_fa)
# check sig2 top model
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_fa, rank = 2,'Fall')
# this does not look as good. stick with nodist weighting.
# but how does sig5 look?
buildBestLandcoverMod(mod_table = coccinellidae_mod_table_fa, rank = 10,'Fall')
# completely different. Must be overfitting.
# reminder - fall coccinellidae densities are mostly zero:
pred_data %>%
  filter(Season == 'Fall', Taxa == 'Coccinellidae') %>%
  ggplot(aes(x = LogDensity, y = Site)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Fall aphids, Log+1 Transformation') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 6)

# maybe could do some zero-inf binomial models, but for now:

# make global mod table (single distweight)
coccinellidae_global_fa <- makeGlobalLandcoverModTab(coccinellidae_mod_table_fa, 1)
# plot var importance for global mod table
plotGlobalVarImportance(coccinellidae_global_fa)
# view table
# View(coccinellidae_global_fa)
# reduced dataset leads to vastly different top models.
# build "best" model
coccinellidae_fin_mod_fa <- buildFinalLandcoverMod(coccinellidae_mod_table_fa, coccinellidae_global_fa, 6)
summary(coccinellidae_fin_mod_fa)
plot(allEffects(coccinellidae_fin_mod_fa, residuals = TRUE), main = 'Coccinellidae, final spring model')

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
geocoris_mod_table_fa <- buildLandcoverModTab('Geocoris', fD_fall)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = geocoris_mod_table_fa, rank = 9, 'Fall')
# Top 8 models are null. 9th doesn't look very good.
# Plot var importance charts
plotVarImportance(mod_table = geocoris_mod_table_fa)
# make global mod table (single distweight)
# try with null model
geocoris_global_fa <- makeGlobalLandcoverModTab(geocoris_mod_table_fa, 1)
# View(geocoris_global_fa)# top mod is null
buildBestLandcoverMod(geocoris_global_fa, 2,'Fall')
geocoris_shan_mod_fa <- buildFinalLandcoverMod(geocoris_global_fa, geocoris_global_fa, 2)
summary(geocoris_shan_mod_fa)
# try with landcover vars
geocoris_global_fa <- makeGlobalLandcoverModTab(geocoris_mod_table_fa, 9)
buildBestLandcoverMod(geocoris_global_fa, 2)
# geocoris_fin_mod_fa <- buildFinalLandcoverMod(geocoris_mod_table_fa, geocoris_global_fa, 2)
# NOT FITTING
# summary(geocoris_fin_mod_fa)
# plot(allEffects(geocoris_fin_mod_fa, residuals = TRUE), main = 'Geocoris, final spring model')
# best mod seems to be shan only


### Aphids ~ landcover ####
#### Spring ####
##### Acyrthosiphon ####

# Build global model selection table
acyrthosiphon_mod_table_sp <- buildLandcoverModTab('Acyrthosiphon', fD_spring)
# Plot best mod and call summary
buildBestLandcoverMod(mod_table = acyrthosiphon_mod_table_sp, rank = 1, 'Spring')
# similar effect with sig1
buildBestLandcoverMod(mod_table = acyrthosiphon_mod_table_sp, rank = 2, 'Spring')
# Plot var importance charts
plotVarImportance(mod_table = acyrthosiphon_mod_table_sp)
# make global mod table (single distweight)
acyrthosiphon_global_sp <- makeGlobalLandcoverModTab(acyrthosiphon_mod_table_sp, 1)
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
# View(acyrthosiphon_fin_mod_sp)

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
  cbind(., to.dummy(.$Treatment, 'Trt')) %>%
  select(2:77)


names(sem_data)


detach("package:lmerTest", unload=TRUE) # this fucks with psem for some reason
spring_acy_sem <- psem(
  lmer(Coccinellidae ~ `weedyWet_sig4` + `Trt.Sham` + (1|site),
     data = sem_data),
  lmer(Acyrthosiphon ~ Coccinellidae + `Trt.Sham`+ (1|site),
     data = sem_data)
  )
plot(spring_acy_sem)
plot(spring_acy_sem, show = 'unstd') #what do standardized coefs mean?


