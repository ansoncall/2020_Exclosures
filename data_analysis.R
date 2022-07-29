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
# Collapse classes
landcover %<>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedyWet = class5 + class8,
         water = class11
         ) %>%
  select(-(class0:class11))

# lengthen
landcover_long <- landcover %>%
  pivot_longer(alfalfa:water,
               names_to = 'class',
               values_to = 'areaScore') %>%
  select(siteId, distanceWeight, site, field, class, areaScore)

ggplot(landcover_long, aes(x = class, y = areaScore, fill = class)) +
  geom_bar( stat = "identity") +
  facet_wrap(~site*distanceWeight, scales = 'free')


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
# ## Aphids ~ predators ####
# ### spring data only ####
#
# # define list of aphid taxa
# short_aphlist <- c('Acyrthosiphon',
#              'Aphis',
#              'Therioaphis')
#
# # create empty lists
# dredges <- list()
#
# for (i in 1:length(short_aphlist)) {
#
#   aph_density <- as.name(short_aphlist[[1]])
#   mGlobal <- mean_density %>%
#     filter(Season=='Spring') %$%
#     lmer(log(eval(aph_density) + 1) ~ (log(Arachnida + 1) +
#                                          log(Coccinellidae + 1) +
#                                          log(Ichneumonidae + 1) +
#                                          log(Nabis + 1) +
#                                          log(Geocoris + 1) +
#                                          log(Anthocoridae + 1) +
#                                          (1|Site) +
#                                          (1|Field)),
#          na.action = 'na.fail', # !!!WATCH OUT ####
#          # UNEXPECTED BEHAVIOR - models for acyrthosiphon fit to subset of data
#          # where all other values are also NA. Effectively throwing away
#          # good data!!
#          REML = FALSE)
#   dredges[[1]] <- dredge(mGlobal) # dredge and store output in list
#
# }
# names(dredges) <- short_aphlist
#
# importance_tab <- lapply(dredges, function (x) {
#   sw(x) %>%
#     tibble(names = names(.),
#            .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
#     mutate(ExpVars = fct_recode(names,
#                                 Anthocoridae = 'log(Anthocoridae + 1)',
#                                 Arachnida = 'log(Arachnida + 1)',
#                                 Coccinellidae = 'log(Coccinellidae + 1)',
#                                 Geocoris = 'log(Geocoris + 1)',
#                                 Ichneumonoidea = 'log(Ichneumonidae + 1)',
#                                 Nabis = 'log(Nabis + 1)'),
#            VarWeights = sw,
#            .keep = 'none') %>%
#     arrange(ExpVars)
# })
#
# names(importance_tab) <- short_aphlist
# # make tables of model selection results
# springTabs <- lapply(dredges, slice_head, n = 5)
# names(springTabs) <- short_aphlist
# # show tables
# springTabs
#
# # extract and examine the top Acyrthosiphon model.
# best.acy.mod <- get.models(springTabs[[1]], subset = 1)[[1]]
# summary(best.acy.mod)
# # must remake to plot effects.
# best.acy.mod <- lmer(log(Acyrthosiphon + 1) ~ log(Coccinellidae + 1) +
#                                                      (1|Site) +
#                                                      (1|Field),
#                      data = mean_density %>% filter(Season == 'Spring'),
#                      REML = FALSE, na.action = 'na.fail')
# summary(best.acy.mod)
# # Note: these two mods are not equal, thanks to 'na.fail' option in for loop.
# # See !!!WATCH OUT note above (line 244)
#
# # make plot
# # png('spring_acy_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Spring')
# # dev.off()
#
# # extract and examine the top Aphis model.
# best.aphis.mod <- get.models(springTabs[[2]], subset = 1)[[1]]
# summary(best.aphis.mod)
# # Must remake to plot effects.
# best.aphis.mod <- lmer(log(Aphis + 1) ~ log(Ichneumonidae + 1) +
#                          (1|Site) +
#                          (1|Field),
#                        data = mean_density %>% filter(Season == 'Spring'),
#                        REML = FALSE)
# # make plot
# # png('spring_aphis_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# plot(allEffects(best.aphis.mod, residuals = TRUE),
#      main = 'Aphis, Spring')
# # dev.off()
#
# # extract and examine the top Therioaphis model.
# best.therio.mod <- get.models(springTabs[[3]], subset = 1)[[1]]
# summary(best.therio.mod)
# # Must remake to plot effects.
# best.therio.mod <- lmer(log(Therioaphis + 1) ~ log(Geocoris + 1) +
#                           (1|Site) +
#                           (1|Field),
#                         data = mean_density %>% filter(Season == 'Spring'),
#                         REML = FALSE)
# # make plot
# # png('spring_therio_effect.jpg',
# #     width = 7,
# #     height = 5,
# #     units = 'in',
# #     res = 300)
# plot(allEffects(best.therio.mod, residuals = TRUE),
#      main = 'Therioaphis, Spring')
# # dev.off()
#
# # make variable importance heatmap
# p <- bind_rows(importance_tab, .id = 'Taxon') %>%
#   mutate(VarWeights = as.numeric(VarWeights),
#          ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
#          Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
#   filter(ExpVars %in% c("Coccinellidae", "Geocoris", "Ichneumonoidea")) %>%
#   ggplot(aes(x = ExpVars, y = Taxon, fill = VarWeights)) +
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0),
#         legend.title = element_text('Variable Importance')) +
#   scale_fill_gradient(low="blue", high="red") +
#   labs(x = '\"Beneficial\" taxon', y = 'Aphid genus', title = 'Spring') +
#   theme(text = element_text(size = 15),
#         axis.text.x=element_text(angle=30,hjust=0.9))
# p
# # ggsave('spring_aphid_varweights.jpg', width = 6, height = 4)
#
# ### fall data only ####
#
# # create empty lists
# dredges <- list()
#
# for (i in 1:length(short_aphlist)) {
#
#   aph_density <- as.name(short_aphlist[[i]])
#   mGlobal <- mean_density %>%
#     filter(Season=='Fall') %$%
#     lmer(log(eval(aph_density) + 1) ~ (log(Arachnida + 1) +
#                                          log(Coccinellidae + 1) +
#                                          log(Ichneumonidae + 1) +
#                                          log(Nabis + 1) +
#                                          log(Geocoris + 1) +
#                                          log(Anthocoridae + 1) +
#                                          (1|Site) +
#                                          (1|Field)),
#          na.action = 'na.fail',
#          REML = FALSE)
#   dredges[[i]] <- dredge(mGlobal) # dredge and store output in list
#
# }
# names(dredges) <- short_aphlist
#
# importance_tab <- lapply(dredges, function (x) {
#   sw(x) %>%
#     tibble(names = names(.),
#            .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
#     mutate(ExpVars = fct_recode(names,
#                                 Anthocoridae = 'log(Anthocoridae + 1)',
#                                 Arachnida = 'log(Arachnida + 1)',
#                                 Coccinellidae = 'log(Coccinellidae + 1)',
#                                 Geocoris = 'log(Geocoris + 1)',
#                                 Ichneumonoidea = 'log(Ichneumonidae + 1)',
#                                 Nabis = 'log(Nabis + 1)'),
#            VarWeights = sw,
#            .keep = 'none') %>%
#     arrange(ExpVars)
# })
# names(importance_tab) <- short_aphlist
# # make tables of model selection results
# fallTabs <- lapply(dredges, slice_head, n = 5)
# names(fallTabs) <- short_aphlist
# # show tables
# # fallTabs
#
# # Extract and examine the top Acyrthosiphon model.
# best.acy.mod <- get.models(fallTabs[[1]], subset = 1)[[1]]
# summary(best.acy.mod)
# # Must remake to plot effects.
# best.acy.mod <- lmer(log(Acyrthosiphon + 1) ~ log(Anthocoridae + 1) +
#                        log(Ichneumonidae + 1) +
#                        (1|Site) +
#                        (1|Field),
#                      data = mean_density %>% filter(Season == 'Spring'),
#                      REML = FALSE)
# # make plot
# # png('fall_acy_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
# plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Fall')
# # dev.off()
#
# # Extract and examine the top Aphis model.
# best.aphis.mod <- get.models(fallTabs[[2]], subset = 1)[[1]]
# summary(best.aphis.mod)
# # Must remake to plot effects.
# best.aphis.mod <- lmer(log(Aphis + 1) ~ log(Arachnida + 1) +
#                        log(Ichneumonidae + 1) +
#                        (1|Site) +
#                        (1|Field),
#                      data = mean_density %>% filter(Season == 'Spring'),
#                      REML = FALSE)
# # make plot
# # png('fall_aphis_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
# plot(allEffects(best.aphis.mod, residuals = TRUE), main = 'Aphis, Fall')
# # dev.off()
#
# # Extract and examine the top Therioaphis model.
# best.therio.mod <- get.models(fallTabs[[3]], subset = 1)[[1]]
# summary(best.therio.mod)
# # Must remake to plot effects.
# best.therio.mod <- lmer(log(Therioaphis + 1) ~ log(Arachnida + 1) +
#                        (1|Site) +
#                        (1|Field),
#                      data = mean_density %>% filter(Season == 'Spring'),
#                      REML = FALSE)
# # make plot
# # png('fall_therio_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
# plot(allEffects(best.therio.mod,
#                 residuals = TRUE),
#      main = 'Therioaphis, Fall')
# # dev.off()
#
# # make plot
# p <- bind_rows(importance_tab, .id = 'Taxon') %>%
#   mutate(VarWeights = as.numeric(VarWeights),
#          ExpVars = fct_reorder(ExpVars, VarWeights, mean, .desc = TRUE),
#          Taxon = fct_reorder(Taxon, .x = VarWeights, .fun = mean)) %>%
#   filter(ExpVars %in% c("Arachnida", "Geocoris", "Ichneumonoidea")) %>%
#   ggplot(aes(x = ExpVars, y = Taxon, fill = VarWeights)) +
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0),
#         legend.title = element_text('Variable Importance')) +
#   scale_fill_gradient(low="blue", high="red") +
#   labs(x = '\"Beneficial\" taxon', y = 'Aphid genus', title = 'Fall') +
#   theme(text = element_text(size = 15),
#         axis.text.x=element_text(angle=30,hjust=0.9))
#
# p
# # ggsave('fall_aphid_varweights.jpg', width = 6, height = 4)
#
# # clean environment
# rm(best.acy.mod, best.aphis.mod, best.therio.mod, dredges, fallTabs,
#    importance_tab, mGlobal, springTabs, aph_density, i)

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

### spring data only ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')

#### formulas ####
# formula for calculating model selection table for predators
buildPredModTab <- function(taxon = 'empty', data = 'empty'){

  if (taxon == 'empty' | !is_tibble(data)) {
    cat(red("Please specify taxon and data \n"))
    stop()
  } else {

    cat(yellow('Taxon:'),
        green(taxon),
        yellow('Data:'),
        green(deparse(substitute(data))),
        '\n')

    mod0 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_no + naturalArid_no + dirt_no + ag_no',
             '+ impermeable_no + weedyWet_no + water_no + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod1 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_const + naturalArid_const + dirt_const + ag_const',
             '+ impermeable_const + weedyWet_const + water_const + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod2 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1',
             '+ impermeable_sig1 + weedyWet_sig1 + water_sig1 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod3 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2',
             '+ impermeable_sig2 + weedyWet_sig2 + water_sig2 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod4 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2',
             '+ impermeable_sig2 + weedyWet_sig2 + water_sig2 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod5 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3',
             '+ impermeable_sig3 + weedyWet_sig3 + water_sig3 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod6 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4',
             '+ impermeable_sig4 + weedyWet_sig4 + water_sig4 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')

    mod7 <- lmer(as.formula(
      paste0('log(',
             taxon,
             ' + 1) ~ alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5',
             '+ impermeable_sig5 + weedyWet_sig5 + water_sig5 + (1|Site)'
      )),
      data = data,
      REML = FALSE,
      na.action = 'na.fail')
    # list global models
    cand_mods <- list(mod0,
                      mod1,
                      mod2,
                      mod3,
                      mod4,
                      mod5,
                      mod6,
                      mod7)
    dfname <- as.name(deparse(substitute(data)))
    for (i in 1:length(cand_mods)) {

      cand_mods[[i]]@call$data <- dfname

    }
    # dredge all models in the list

    dredges <- lapply(cand_mods, dredge, )
    # Rbind the elements of the list together. This forces recalculation of AICc
    mod_table <- rbind(dredges[[1]],
                       dredges[[2]],
                       dredges[[3]],
                       dredges[[4]],
                       dredges[[5]],
                       dredges[[6]],
                       dredges[[7]],
                       dredges[[8]])
  }
}

# extract and examine models by rank
plotBestPredMod <- function(mod_table, rank) {
  best.mod <- get.models(mod_table, subset = rank)[[1]]
  getPred <- best.mod@call$formula[[2]]
  cat(green('Summary for model',  'rank =', rank, '\n'))
  cat(yellow('Warning: more than one plot may be printed \n'))
  plot(allEffects(best.mod, residuals = TRUE),
       main = paste(paste(getPred, collapse = ''),
                     'landcover model rank',
                     as.character(rank)))
  summary(best.mod)
}

# plot the variable importance
plotPredVarImportance <- function (mod_table){

  getPred <- get.models(mod_table, subset = 1)[[1]]@call$formula[[2]]
  importance_tab <- sw(mod_table) %>% #tibble(names = names(.))
    tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
    arrange(names) %>%
    separate(names, c('class', 'distWeight'), sep = "_") %>%
    mutate(distWeight = as_factor(recode(distWeight,
                                         `no` = 'no decay',
                                         `const` = 'constant',
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
         title = paste(paste(getPred, collapse = ''), 'variable importance'))

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
makeGlobalPredModTab <- function(mod_table, rank = 1){
  mod <- get.models(mod_table, subset = rank)[[1]]
  # Global model is rank-deficient. We will have to 'trick' dredge.

  # mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
  #                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
  #                     (1|Site),
  #                  data = fD_spring_sub) # This doesn't work.
  response <- names(mod@frame) %>% head(1)
  predictors <- names(mod@frame) %>% head(-1) %>% tail(-1)
  vars.all <- c(predictors,
                'shan',
                'rich',
                'total_cover')
  # Write formula for full global model.
  form <- formula(paste0(response, '~',
                         paste0(vars.all, collapse='+'),
                         '+(1|Site)'))

  # Prepare formula with reduced number of terms.
  vars.red <- c('shan', 'rich', 'total_cover')
  form.red <- formula(paste0('log(Arachnida + 1) ~',
                             paste0(vars.red,collapse='+'),
                             '+(1|Site)'))
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
plotPredGlobalVarImportance <- function(global_tab){
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
buildFinPredMod <- function(no_margin_tab, global_tab, rank = 1){

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



#### Arachnida ####
# Build global model selection table
arachnida_mod_table <- buildPredModTab('Arachnida', fD_spring)
# Plot best mod and call summary
plotBestPredMod(mod_table = arachnida_mod_table, rank = 1)
# Plot var importance charts
plotPredVarImportance(mod_table = arachnida_mod_table)
# make global mod table (single distweight)
arachnida_global <- makeGlobalPredModTab(arachnida_mod_table, 1)
# plot var importance for global mod table
plotPredGlobalVarImportance(arachnida_global)
# build best model
arachnida_fin_mod <- buildFinPredMod(arachnida_mod_table, arachnida_global, 1) # Error
summary(arachnida_fin_mod)
plot(allEffects(arachnida_fin_mod, residuals = TRUE), main = 'Arachnida, final model')


#### Coccinellidae ####
# Build global model selection table
coccinellidae_mod_table <- buildPredModTab('Coccinellidae', fD_spring)
# Plot best mod and call summary
plotBestPredMod(mod_table = coccinellidae_mod_table, rank = 1)
# Plot var importance charts
plotPredVarImportance(mod_table = coccinellidae_mod_table)
# make global mod table (single distweight)
coccinellidae_global <- makeGlobalPredModTab(coccinellidae_mod_table, 1)
# plot var importance for global mod table
plotPredGlobalVarImportance(coccinellidae_global)
# build best model
coccinellidae_fin_mod <- buildFinPredMod(coccinellidae_mod_table, coccinellidae_global, 1) # Error
summary(coccinellidae_fin_mod)
plot(allEffects(coccinellidae_fin_mod, residuals = TRUE), main = 'Coccinellidae, final model')

##### Compare old and new ####
# 2022-07-28T17:39:44Z this is old now bc of superdove data. Deleted.

#### Ichneumonidae ####
# Build global model selection table
ichneumonidae_mod_table <- buildPredModTab('Ichneumonidae', fD_spring)
# Plot best mod and call summary
plotBestPredMod(mod_table = ichneumonidae_mod_table, rank = 1)
# Plot var importance charts
plotPredVarImportance(mod_table = ichneumonidae_mod_table)
# make global mod table (single distweight)
ichneumonidae_global <- makeGlobalPredModTab(ichneumonidae_mod_table, 1)
# plot var importance for global mod table
plotPredGlobalVarImportance(ichneumonidae_global)
# build best model
ichneumonidae_fin_mod <- buildFinPredMod(ichneumonidae_mod_table,
                                         ichneumonidae_global,
                                         2)
summary(ichneumonidae_fin_mod)
plot(allEffects(ichneumonidae_fin_mod, residuals = TRUE), main = 'Ichneumonidae, final model')
# hypothetical ich model with dirt_const and shan:
ich_test_mod <- lmer(log(Ichneumonidae + 1) ~ dirt_const + shan + (1|Site),
                     data = fD_spring_sub,
                     REML = FALSE)
summary(ich_test_mod)
# this is model 6 in the 'global' mod table btw
ich_from_tab <-get.models(ichneumonidae_global, subset = 6)[[1]]
summary(ich_from_tab)

# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_no +
               naturalArid_no +
               dirt_no +
               ag_no +
               impermeable_no +
               weedyWet_no +
               water_no +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_const +
               naturalArid_const +
               dirt_const +
               ag_const +
               impermeable_const +
               weedyWet_const +
               water_const +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_sig1 +
               naturalArid_sig1 +
               dirt_sig1 +
               ag_sig1 +
               impermeable_sig1 +
               weedyWet_sig1 +
               water_sig1 +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_sig2 +
               naturalArid_sig2 +
               dirt_sig2 +
               ag_sig2 +
               impermeable_sig2 +
               weedyWet_sig2 +
               water_sig2 +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_sig3 +
               naturalArid_sig3 +
               dirt_sig3 +
               ag_sig3 +
               impermeable_sig3 +
               weedyWet_sig3 +
               water_sig3 +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_sig4 +
               naturalArid_sig4 +
               dirt_sig4 +
               ag_sig4 +
               impermeable_sig4 +
               weedyWet_sig4 +
               water_sig4 +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod6 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa_sig5 +
               naturalArid_sig5 +
               dirt_sig5 +
               ag_sig5 +
               impermeable_sig5 +
               weedyWet_sig5 +
               water_sig5 +
               (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5,
                  mod6)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]],
                   dredges[[7]])
# Print the table.
# View(mod_table)

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)

plot(allEffects(best.mod, residuals = TRUE),
     main = 'Ichneumonidae, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = '_') %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Ichneumonidae + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c(
  'alfalfa_const',
  'naturalArid_const',
  'dirt_const',
  'ag_const',
  'impermeable_const',
  'weedyWet_const',
  'water_const',
  'shan',
  'rich',
  'total_cover')
# Write formula for full global model.
form <- formula(paste0('log(Ichneumonidae + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Ichneumonidae + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 4), eval))
summary(ms1)
best.mod <- get.models(ms1, subset = 2)[[1]]
summary(best.mod)
# best mod here is null
# second best mod is 'impermeable_const' p=0.0375
# if using 'sig4' weighting, the second best mod is 'shan' p=0.0496
# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Ichneumonidae + 1) ~
                  impermeable_const + (1|Site),
                data = fD_spring_sub)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Ichneumonidae, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Geocoris ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Geocoris + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Geocoris + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Geocoris + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Geocoris + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Geocoris + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Geocoris + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
# View(mod_table)

# Extract and examine the best model.
# Top 6 mods are null. 7th best is shown below.
best.mod <- get.models(mod_table, subset = 7)[[1]]
summary(best.mod)
plot(allEffects(best.mod, residuals = TRUE),
     main = 'Geocoris, best landcover model')


png('spring_geo_effect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(allEffects(best.mod, residuals = TRUE),
     main = '', #Coccinellidae - best landcover model
     partial.residual = list(lwd = 0),
     axes = list(x = list(natural2 =
                            list(lab =
                                   list(label =
                                          'Weighted proportion of \"natural\" landcover',
                                        cex = 1.5))),
                 y = list(lab = list(label = 'log(Geocoris density)', cex = 1.5))))
dev.off()

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance %>% filter(distWeight != 'constant'),
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'weight')
p
ggsave('spring_geo_groupvarweights.jpg', width = 7, height = 5)

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Geocoris + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa2',
              'bare2', 'disturbed2', 'natural2', 'wet2')
# Write formula for full global model.
form <- formula(paste0('log(Geocoris + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Geocoris + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Geocoris + 1) ~ natural2 + total_cover +
                  (1|Site),
                data = fD_spring, REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Geocoris, final model')

# png('spring_geo_natEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('natural2', fin.mod, residuals = TRUE),
     main = NULL,
     xlab = 'Weighted proportion of \"natural\" landcover',
     ylab = 'log(Geocoris density)')
# dev.off()
# png('spring_geo_coverEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('total_cover', fin.mod, residuals = TRUE),
     main = '', #Coccinellidae - best landcover model
     partial.residual = list(lwd = 0),
     axes = list(x = list(total_cover =
                            list(lab =
                                   list(label =
                                          'Mean % cover in field margins',
                                        cex = 1.5))),
                 y = list(lab = list(label = 'log(Geocoris density)', cex = 1.5))))
# dev.off()
# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

### fall data only ####

# Subset data for modeling with landcover classes.
fD_fall <- field_data_grouped %>%
  filter(Season == 'Fall') %>%
  mutate_at(72:101, scale) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')

#### Arachnida ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Arachnida + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Arachnida + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Arachnida + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Arachnida + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Arachnida + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Arachnida + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
# Best model is null.
# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Arachnida, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_fall_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa0',
              'bare0', 'disturbed0', 'natural0', 'wet0')
# Write formula for full global model.
form <- formula(paste0('log(Arachnida + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Arachnida + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
# Note: seems to be overfitting with m.max = 4
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Arachnida + 1) ~
                  wet0 + alfalfa0 + (1|Site),
                data = fD_fall_sub, REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Arachnida, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)




#### Coccinellidae ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Coccinellidae + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Coccinellidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Coccinellidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Coccinellidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Coccinellidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Coccinellidae + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
# Best mod is null.
# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Coccinellidae, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance %>% filter(distWeight != 'constant'),
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'weight')
p
ggsave('fall_cocc_groupvarweights.jpg', width = 7, height = 5)

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Coccinellidae + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_fall_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa0',
              'bare0', 'disturbed0', 'natural0', 'wet0')
# Write formula for full global model.
form <- formula(paste0('log(Coccinellidae + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Coccinellidae + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
# Note: 3 or 4 terms seems to lead to overfitting.
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Coccinellidae + 1) ~
                  natural5 + (1|Site),
                data = fD_fall_sub,
                REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Coccinellidae, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Ichneumonidae ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Ichneumonidae, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Ichneumonidae + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_fall_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa5',
              'bare5', 'disturbed5', 'natural5', 'wet5')
# Write formula for full global model.
form <- formula(paste0('log(Ichneumonidae + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Ichneumonidae + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Ichneumonidae + 1) ~
                  disturbed5 + (1|Site),
                data = fD_fall,
                REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Ichneumonidae, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Geocoris ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Geocoris + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Geocoris + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Geocoris + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Geocoris + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Geocoris + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Geocoris + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Geocoris, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Geocoris + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_fall_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa4',
              'bare4', 'disturbed4', 'natural4', 'wet4')
# Write formula for full global model.
form <- formula(paste0('log(Geocoris + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Geocoris + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Geocoris + 1) ~
                  natural4 + (1|Site),
                data = fD_fall)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Geocoris, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)


### Aphids ~ landcover ####
#### Spring ####
##### Acyrthosiphon ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)

# Fix na.action
options(na.action = 'na.fail')
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
View(mod_table)

# Extract and examine the best model.
# Best mods are all null. 7th best shown below.
best.mod <- get.models(mod_table, subset = 7)[[1]]
summary(best.mod)

# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Acyrthosiphon, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'weight')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Acyrthosiphon + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa0',
              'bare0', 'disturbed0', 'natural0', 'wet0')
# Write formula for full global model.
form <- formula(paste0('log(Acyrthosiphon + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(Acyrthosiphon + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red, data = fD_spring_sub,
                 REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Acyrthosiphon + 1) ~ rich + natural0 + (1|Site),
                data = fD_spring_sub,
                REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Acyrthosiphon, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Fall ####
##### Non-acyrthosiphon #####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(NonAcy + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(NonAcy + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(NonAcy + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(NonAcy + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(NonAcy + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(NonAcy + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE,
             na.action = 'na.fail')
# List global models.
cand_mods <- list(mod0,
                  mod1,
                  mod2,
                  mod3,
                  mod4,
                  mod5)
# Dredge all models in the list.
dredges <- lapply(cand_mods, dredge)
# Rbind the elements of the list together. This forces recalculation of AICc
mod_table <- rbind(dredges[[1]],
                   dredges[[2]],
                   dredges[[3]],
                   dredges[[4]],
                   dredges[[5]],
                   dredges[[6]])
# Print the table.
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
# Best model is null.
# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'NonAcy, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'sw')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(sw))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'weight')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(NonAcy + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_fall_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa5',
              'bare5', 'disturbed5', 'natural5', 'wet5')
# Write formula for full global model.
form <- formula(paste0('log(NonAcy + 1) ~',
                       paste0(vars.all, collapse='+'),
                       '+(1|Site)'))

# Prepare formula with reduced number of terms.
vars.red <- c('shan', 'rich', 'total_cover')
form.red <- formula(paste0('log(NonAcy + 1) ~',
                           paste0(vars.red,collapse='+'),
                           '+(1|Site)'))
# Fit reduced model.
fmod.red <- lmer(form.red,
                 data = fD_fall_sub,
                 REML = FALSE,
                 na.action = 'na.fail')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = sw)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'sw')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(NonAcy + 1) ~
                  alfalfa5 + (1|Site),
                data = fD_fall, REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'NonAcy, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

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
## Spring ####
### Acyrthosiphon ####
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
  select(2:6, 11:84) %>%
  rename(bare2 = )

names(sem_data)

spring_acy_sem <- psem(
  lm(Coccinellidae ~ `2_classification_sig1` +
       `6_classification_sig1` + `Trt.Sham`,
     data = sem_data),
  lm(Acyrthosiphon ~ Coccinellidae + `Trt.Sham`,
     data = sem_data)
  )
