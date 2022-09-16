# Load packages ####
library(car) # for Anova() on lmer model objects
library(corrr) # for correlation plots of landcover vars
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
## CHECK #####
## make sure you are using the correct landcover table!
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
  mutate(LogDensity = log(Density+1)) %>%
  filter(Taxa %in% c('Acyrthosiphon', 'NonAcy'))

# build density plot
ggplot(data = aph_data, aes(y = Taxa, x = LogDensity, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Aphids, Log+1 Transformation',
       x = 'log(Density)',
       y = 'Taxon') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 10)
ggsave('aph-abund.png', width = 12.5, height = 6, units = 'in')
# barchart for slideshow
# spring
ggplot(data = aph_data %>%
         filter(Season == 'Spring',
                Taxa %in% c('Acyrthosiphon', 'NonAcy')),
       aes(x = LogDensity, y = Taxa, fill = Taxa)) +
  geom_boxplot() +
  xlim(c(0, 8.5)) +
  scale_fill_manual(values = c('#548235', '#3b3838'),
                    guide = 'none') +
  labs(y = '', x = 'log(Density)')+
  coord_flip() +
  theme(text = element_text(size = 20)) +
  scale_y_discrete(labels=c("Pea aphids +\n Blue alfalfa aphids","Other taxa"))


ggsave('spring_aphid_density.jpg', width = 7, height = 5)
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
         wet =  class8,
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
# ggsave(filename = 'lc.jpg', width = 9, height = 7, units = 'in')


ggplot(landcover_long %>% filter(site == 'Lovelock'),
       # aes - limit x axis to 4 characters
       aes(x = sub(class, pattern = "(\\w{4}).*", replacement = "\\1."),
           y = areaScore,
           fill = class)) +
  geom_bar( stat = "identity") +
  facet_wrap(ncol = 7, ~field*distanceWeight, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = 'lc.jpg', width = 9, height = 7, units = 'in')

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
             'Therioaphis', 'NonAcy')

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

## better plot
bmod <- lmer(Acyrthosiphon ~ Coccinellidae +
         (1|Site)+(1|Field),
       na.action = 'na.omit',
       REML = FALSE,
       data = mean_density %>%
         filter(Season=='Spring') %>%
         mutate(Coccinellidae = log(Coccinellidae+1),
                Acyrthosiphon = log(Acyrthosiphon+1)))

summary(bmod)

dfb <-mean_density %>%
  filter(Season=='Spring') %>%
  mutate(Coccinellidae = log(Coccinellidae+1),
         Acyrthosiphon = log(Acyrthosiphon+1))
trellis.device()
trellis.par.get()


png('lattice.png', height = 6, width = 6, units = 'in', res = 300)
trellis.par.set(list(par.xlab.text = list(cex=2),
                     par.ylab.text = list(cex=2),
                     par.main.text = list(col = "blue", cex=0.5)))
plot(effect('Coccinellidae',bmod,
                residuals = T),
     partial.residuals = list(smooth=F),
     axes = list(x = list(Coccinellidae = list(lab = 'Log(Ladybug density)')),
                 y = list(lab = 'Log(Aphid density)')),
     main = NULL,
     lattice=list(key.args=list(axis.text = list(cex=10))))

dev.off()
dfb %>%
  ggplot(aes(Coccinellidae, Acyrthosiphon))+
  geom_point()


library(ggeffects)
b <- ggemmeans(bmod, terms = 'Coccinellidae')

ggplot(b, aes(x, predicted)) +
  geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)

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
  # filter(ExpVars %in% c("Coccinellidae", "Geocoris", "Ichneumonoidea")) %>%
  filter(Taxon %in% c("Acyrthosiphon")) %>%
  ggplot(aes(x = reorder(ExpVars, desc(VarWeights)), y = VarWeights, fill = VarWeights)) +
  geom_col() +
  scale_fill_gradient('', low="blue", high="red", breaks = c(0.3,0.9), labels = c('low','high')) +
  labs(x = '', y = 'Variable importance') +
  theme_gray(base_size = 25)+
  theme( legend.position = c(0.95,0.85)) +
  scale_x_discrete(labels=c('Ladybugs', 'Spiders', 'Parasitoid\nwasps', 'Damsel\nbugs', 'Minute pirate\nbugs','Bigeyed\nbugs'))
p
ggsave('spring_aphid_varweights.jpg', width = 12.5, height = 7.5)

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

# Nonacy
# Extract and examine the top NonAcy model.
best.nonacy.mod <- get.models(fallTabs[[4]], subset = 1)[[1]]
summary(best.nonacy.mod)
# Must remake to plot effects.
best.nonacy.mod <- lmer(log(NonAcy + 1) ~ log(Arachnida + 1) +
                          (1|Site) +
                          (1|Field),
                        data = mean_density %>% filter(Season == 'Spring'),
                        REML = FALSE)
# make plot
# png('fall_nonacy_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
plot(allEffects(best.nonacy.mod,
                residuals = TRUE),
     main = 'NonAcy, Fall')

# ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
#   geom_tile() +
#   theme(axis.text.x=element_text(angle = 45, hjust = 0),
#         axis.text.y = element_text(angle = 45)) +
#   scale_fill_gradient(low="blue", high="red") +
#   labs(x = 'Landcover class',
#        y = 'Distance weighting algorithm',
#        title = paste0('log(', getPred[2], ')',
#                       ' Variable importance, ',
#                       season)
#   )

# make plot
p <- bind_rows(importance_tab, .id = 'Taxon') %>%
  mutate(VarWeights = as.numeric(VarWeights),
         ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
         Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
  filter(ExpVars %in% c("Arachnida", "Geocoris", "Ichneumonoidea", "Anthocoridae")) %>%
  ggplot(aes(x = ExpVars, y = Taxon, fill = VarWeights)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        legend.title = element_text('Variable Importance')) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = '\"Beneficial\" taxon', y = 'Aphid genus', title = 'Fall') +
  theme(text = element_text(size = 15),
        axis.text.x=element_text(angle=30,hjust=0.9))

p
ggsave('fall_aphid_varweights.jpg', width = 6, height = 4)

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
    theme_gray(base_size = 20)+
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
    theme_gray(base_size = 20) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_gradient(low="blue", high="red") +
    labs(x = 'Landcover class',
         y = 'Distance weighting algorithm',
         fill = '')

  qq <- ggplotly(q, tooltip = 'weight')
  subplot(pp, qq, widths = c(7/8, 1/8))


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




### spring data only ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')


### Correlation of landcover factors ####
distList <- c('_sig1',
              '_sig2',
              '_sig3',
              '_sig4',
              '_sig5',
              '_const',
              'no')

distList2 <- c('Very Aggressive',
              'Aggressive',
              'Moderately aggressive',
              'Moderate',
              'Slight',
              'Constant',
              'None')

dVarList <- c()
library(tidyselect)
for (i in 1:length(distList)) {

  dVarList[[i]] <- fD_spring %>%
    select(ends_with(distList[[i]])) %>%
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


}

dVarList
for (i in 1:length(distList)) {

  print(rplot(dVarList[[i]], print_cor = TRUE, legend = FALSE, .order = 'alphabet') +
          theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
          labs(title = paste(distList2[[i]]))
  )
  ggsave(paste0('corr', distList[[i]], '.png'),
         width = 4,
         height = 4,
         units = 'in')

}

fD_spring %>%
  select(ends_with('sig3'), id) %>%
  rename_with(~ str_remove(.x, "_.+")) %>%
  ggplot(aes(wet, dirt)) +
  geom_point() +
  geom_text(aes(label = id), hjust = 0, vjust = 0)+
  labs(x = 'Riparian, moderate weighting',
       y = 'Bare soil + dirt road,\nmoderate weighting')+
  coord_cartesian(xlim = c(-1.2, 2.4))

ggsave('dirtByWet_sig3.png', width = 7, height = 5, units = 'in')

### Flooded vs sprinklers ####

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

fD_spring_flood <- cbind(fD_spring, wateringMethod)

floodMod <- lmer(log(Acyrthosiphon+1)~wateringMethod+(1|site),
                 data = fD_spring_flood)
floodMod
summary(floodMod)

plot(allEffects(floodMod))

floodMod2 <- lm(log(Acyrthosiphon+1)~wateringMethod, data = fD_spring_flood)
floodMod2
summary(floodMod2)

plot(allEffects(floodMod2))

p=ggplot(fD_spring_flood,
       aes(x = wateringMethod,
           y = log(Acyrthosiphon + 1))) +
  # geom_boxplot() +
  # geom_jitter() +
  # geom_violin()
  geom_jitter(aes(color = Site), size = 6, width = 0.06) +
  # geom_smooth(method = 'lm') +
  stat_summary(fun.data = 'mean_cl_boot') +
  labs(title = 'mean_cl_boot',
       subtitle = 'basic nonparametric bootstrap for obtaining confidence limits')

q=ggplot(fD_spring_flood,
       aes(x = wateringMethod,
           y = log(Acyrthosiphon + 1))) +
  # geom_boxplot() +
  # geom_jitter() +
  # geom_violin()
  geom_jitter(aes(color = Site), size = 6, width = 0.06) +
  # geom_smooth(method = 'lm') +
  stat_summary(fun.data = 'mean_cl_normal') +
  labs(title = 'mean_cl_normal',
       subtitle = 'lower and upper Gaussian confidence limits based on the t-distribution')

r=ggplot(fD_spring_flood,
       aes(x = wateringMethod,
           y = log(Acyrthosiphon + 1))) +
  geom_boxplot() +
  # geom_jitter() +
  # geom_violin()
  geom_jitter(aes(color = Site), size = 6, width = 0.06) +
  # geom_smooth(method = 'lm') +
  # stat_summary(fun.data = 'mean_cl_normal') +
  labs(title = 'boxplot',
       subtitle = 'a normal boxplot')

grid.arrange(p,q,r, nrow =1)


flooded <- fD_spring_flood %>% filter(wateringMethod=='Flooded') %>%
  mutate(Acyrthosiphon = log(Acyrthosiphon+1)) %>%
  select(Acyrthosiphon)
sprinklers <- fD_spring_flood %>% filter(wateringMethod=='Sprinklers') %>%
  mutate(Acyrthosiphon = log(Acyrthosiphon+1)) %>%
  select(Acyrthosiphon)

t.test(flooded, sprinklers, var.equal = TRUE)
t.test(flooded, sprinklers)

final_fig <- ggplot(fD_spring_flood,
                    aes(x = wateringMethod,
                        y = log(Acyrthosiphon + 1))) +
  geom_jitter(aes(color = Site), size = 6, width = 0.1) +
  stat_summary(fun.data = 'mean_cl_boot', geom="errorbar", width = 0.3)+
  labs(
       y = 'log(Acyrthosiphon density)',
       x = 'Irrigation method')+
  theme_classic(base_size = 20) +
  scale_color_brewer(type = 'qual', palette = 2)
final_fig
ggsave('watering.png', width = 5.5, height = 6.5)



### sum of nat ####
sums <- field_data %>%
  filter(Season == 'Spring') %>% # must not use scaled data
  mutate(sum_no = rowSums(across(ends_with('no')))) %>%
  mutate(sum_sig4 = rowSums(across(ends_with('sig4')))) %>%
  mutate(sum_sig2 = rowSums(across(ends_with('sig2')))) %>%
  select(id, ends_with('no'), ends_with('sig4'), ends_with('sig2')) %>%
  select(id, starts_with('naturalArid'), starts_with('sum')) %>%
  mutate(ratio_no = (naturalArid_no/sum_no)*100,
         ratio_sig4 = (naturalArid_sig4/sum_sig4)*100,
         ratio_sig2 = (naturalArid_sig2/sum_sig2)*100) %>%
  select(id, starts_with('ratio')) %>%
  filter(id == 'Minden02')


sums
ratios<-sums %>% select(-id)
percentile <- ratios$ratio_no
ggplot() +
  geom_col(aes("", 100)) +
  geom_col(aes("", percentile), fill = "#CF3FFF") +
  coord_flip() +
  theme_void()
ggsave('naNo.png', width = 3, height = 0.5, units = 'in')
percentile <- ratios$ratio_sig4
ggplot() +
  geom_col(aes("", 100)) +
  geom_col(aes("", percentile), fill = "#CF3FFF") +
  coord_flip() +
  theme_void()
ggsave('sig4No.png', width = 3, height = 0.5, units = 'in')
percentile <- ratios$ratio_sig2
ggplot() +
  geom_col(aes("", 100)) +
  geom_col(aes("", percentile), fill = "#CF3FFF") +
  coord_flip() +
  theme_void()
ggsave('sig2No.png', width = 3, height = 0.5, units = 'in')


#### General figs ####
# define list of taxa of interest
taxaList <- c('Arachnida',
           'Coccinellidae',
           'Geocoris',
           'Ichneumonidae',
           'Acyrthosiphon')

# read in tabs
tabList <- c()
for (i in 1:length(taxaList)) {

  tabList[[i]] <- readRDS(paste0('modTabs/', taxaList[[i]], '_spring_8class'))

}
# generate vi heatmaps ####
viList <- c()
for (i in 1:length(taxaList)) {

  viList[[i]] <- plotVarImportance(tabList[[i]], 'Spring')

}
viList[[5]]

ggsave('cocc_var_heatmap.png', width = 12.5, height = 7, units = 'in')

# view head of tabs
heads <- c()
for (i in 1:length(taxaList)) {

  heads[[i]] <- head(tabList[[i]])

}

viList[[5]]
View(tabList[[5]])

# plot top mod effects (one at a time, manually)

k =5
  plotName <- paste0(taxaList[[k]], '_rank1_mod_effect.png')
  plotMod <- get.models(tabList[[k]], 2)[[1]]
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
library(DiagrammeRsvg)
plot(spring_acy_sem)%>%
  export_svg %>% charToRaw %>% rsvg_png("graph.png")
install.packages('rsvg')
library(rsvg)
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




