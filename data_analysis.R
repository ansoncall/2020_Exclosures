# Load packages ####
library(broom.mixed) # for augment() and tidy() to easily build tibbles
library(car) # for Anova() on lmer model objects
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
landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffddddddddc')
vegPlots <- read_csv('tidy_data/vegPlots.csv', col_types = 'fffffff')
vegSites <- read_csv('tidy_data/vegSites.csv', col_types = 'f')

# define convenience identifiers ####
# define lists of predator and aphid taxa
predlist <- c('Arachnida','Coccinellidae','Ichneumonidae',
              'Nabis', 'Geocoris', 'Anthocoridae')
aphlist <- c('Acyrthosiphon', 'Aphis', 'Therioaphis', 'AllAph', 'NonAcy')

# define predator- and aphid-only tibbles
# filter data to include only predators
preddata <- right_join(data_long %>% filter(Treatment != 'Pre-'),
                       data_long %>% filter(Treatment == 'Pre-'),
                       by = c('Site','Field','Plot','Taxa','Season')) %>%
  select(Site:Plot,
         Treatment = Treatment.x,
         Taxa,
         Density = Density.x,
         Pre_Density = Density.y,
         Season) %>%
  filter(Taxa %in% predlist)

aphdata <- right_join(data_long %>% filter(Treatment != 'Pre-'),
                      data_long %>% filter(Treatment == 'Pre-'),
                      by = c('Site','Field','Plot','Taxa','Season')) %>%
  select(Site:Plot,
         Treatment = Treatment.x,
         Taxa,
         Density = Density.x,
         Pre_Density = Density.y,
         Season) %>%
  filter(Taxa %in% aphlist)

# make df of total abundance across treatments
mean_density <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  separate(id, c('Site', 'Field', 'Plot', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density)


# Aphid histogram ####
aph_data <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  filter(Taxa %in% aphlist) %>%
  mutate(LogDensity = log(Density+1))

# build density plot app
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
ggsave('fall_aphid_density.jpg', width = 7, height = 5)

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
ggsave('total_aphid_density.jpg', width = 1.5, height = 4.5)

## Notes ####
## Just make this a boxplot - they are more familiar. Also consider
## putting season data side by side and faceting on taxon. Also need consistent
## colors.

# Predator histogram ####
pred_data <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  filter(Taxa %in% predlist) %>%
  mutate(LogDensity = log(Density+1), Taxa = fct_relevel(Taxa, sort))

# build density plot app
ggplot(data = pred_data, aes(y = Taxa, x = LogDensity, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = 'Predators, Log+1 Transformation') +
  theme(legend.position = 'none') +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 6)
# barchart for slideshow
ggplot(data = pred_data,
       aes(x = LogDensity, y = Taxa, fill = Season)) +
  geom_boxplot() +
  coord_flip() +
  xlim(c(0, 6)) +
  labs(y = 'Predator taxon', x = 'log(Density)') +
  scale_fill_manual(values = c('#008000','#993300')) +
  theme(text = element_text(size = 20),
        axis.text.x=element_text(angle=30,hjust=0.9))
ggsave('spring_pred_density.jpg', width = 10, height = 5)
ggplot(data = pred_data %>%
         filter(Season == 'Fall'),
       aes(x = LogDensity, y = Taxa, fill = Taxa)) +
  geom_boxplot() +
  coord_flip() +
  xlim(c(0, 6)) +
  guides(fill = 'none') +
  labs(y = 'Predator taxon', x = 'log(Density)')
ggsave('fall_pred_density.jpg', width = 12, height = 5)

## demo figure for field effect on density
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
ggsave('example_field_effect.jpg', width = 4, height = 5)
## Notes ####
## Just make this a boxplot - they are more familiar. Also consider
## putting season data side by side and faceting on taxon. Also need consistent
## colors.

# Landcover data summary ####
# lengthen
landcover_long <- landcover %>%
  pivot_longer(`0_classification`:`7_classification`,
               names_to = 'class',
               values_to = 'areaScore') %>%
  select(id, distanceWeight, site, field, class, areaScore)

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
## Aphids ~ predators ####
# spring data only

aphlist <- c('Acyrthosiphon',
             'Aphis',
             'Therioaphis')
aphDFs <- list() # create empty lists
dredges <- list()

for (i in 1:length(aphlist)) {

  aph_density <- as.name(aphlist[[i]])
  mGlobal <- mean_density %>%
    filter(Season=='Spring') %$%
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
names(dredges) <- aphlist

importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tidy() %>%
    mutate(ExpVars = fct_recode(names,
                                Anthocoridae = 'log(Anthocoridae + 1)',
                                Arachnida = 'log(Arachnida + 1)',
                                Coccinellidae = 'log(Coccinellidae + 1)',
                                Geocoris = 'log(Geocoris + 1)',
                                Ichneumonoidea = 'log(Ichneumonidae + 1)',
                                Nabis = 'log(Nabis + 1)'),
           VarWeights = x,
           .keep = 'none') %>%
    arrange(ExpVars)
})

names(importance_tab) <- aphlist
# make kable
springTabs <- lapply(dredges, slice_head, n = 5)
names(springTabs) <- aphlist

# Extract and examine the top Acyrthosiphon model.
best.acy.mod <- get.models(springTabs[[1]], subset = 1)[[1]]
summary(best.acy.mod)
# Must remake to plot effects.
best.acy.mod <- lmer(log(Acyrthosiphon + 1) ~ log(Coccinellidae + 1) +
                                                     (1|Site) +
                                                     (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'))
plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Spring')
png('spring_acy_effect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Spring')
dev.off()

# Extract and examine the top Therioaphis model.
best.therio.mod <- get.models(springTabs[[3]], subset = 1)[[1]]
summary(best.therio.mod)
# Must remake to plot effects.
best.therio.mod <- lmer(log(Therioaphis + 1) ~ log(Geocoris + 1) +
                       (1|Site) +
                       (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'))
plot(allEffects(best.therio.mod, residuals = TRUE), main = 'Therioaphis, Spring')
png('spring_therio_effect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(allEffects(best.therio.mod, residuals = TRUE), main = 'Therioaphis, Spring')
dev.off()


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

ggsave('spring_aphid_varweights.jpg', width = 6, height = 4)

# fall data only

aphDFs <- list() # create empty lists
dredges <- list()

for (i in 1:length(aphlist)) {

  aph_density <- as.name(aphlist[[i]])
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
names(dredges) <- aphlist

importance_tab <- lapply(dredges, function (x) {
  sw(x) %>%
    tidy() %>%
    mutate(ExpVars = fct_recode(names,
                                Anthocoridae = 'log(Anthocoridae + 1)',
                                Arachnida = 'log(Arachnida + 1)',
                                Coccinellidae = 'log(Coccinellidae + 1)',
                                Geocoris = 'log(Geocoris + 1)',
                                Ichneumonoidea = 'log(Ichneumonidae + 1)',
                                Nabis = 'log(Nabis + 1)'),
           VarWeights = x,
           .keep = 'none') %>%
    arrange(ExpVars)
})
names(importance_tab) <- aphlist
# make kable
fallTabs <- lapply(dredges, slice_head, n = 5)
names(fallTabs) <- aphlist


# Extract and examine the top Acyrthosiphon model.
best.acy.mod <- get.models(fallTabs[[1]], subset = 1)[[1]]
summary(best.acy.mod)
# Must remake to plot effects.
best.acy.mod <- lmer(log(Acyrthosiphon + 1) ~ log(Anthocoridae + 1) +
                       log(Ichneumonidae + 1) +
                       (1|Site) +
                       (1|Field),
                     data = mean_density %>% filter(Season == 'Spring'))
plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Fall')
png('fall_acy_effect.jpg', width = 14, height = 5, units = 'in', res = 300)
plot(allEffects(best.acy.mod, residuals = TRUE), main = 'Acyrthosiphon, Fall')
dev.off()

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

ggsave('fall_aphid_varweights.jpg', width = 6, height = 4)

## Predators ~ ####
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
              values_from = `0_classification`:`7_classification`) %>%
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

# combine classes
# // key:
# // 0. Alfalfa
# // 1. Disturbed_arid
# // 2. Bare_soil (non-road)
# // 3. Dirt_road
# // 4. Impermeable (paved, concrete, rooftops, etc.)
# // 5. Natural_riparian
# // 6. Natural_arid
# // 7. Wet_ditches

# new key
# alfalfa: 0
# disturbed: 1
# bare: 234
# natural: 5
# wet: 67

field_data_grouped <- field_data %>%
  mutate(alfalfa0 = `0_classification_const`,
         alfalfa1 = `0_classification_sig1`,
         alfalfa2 = `0_classification_sig2`,
         alfalfa3 = `0_classification_sig3`,
         alfalfa4 = `0_classification_sig4`,
         alfalfa5 = `0_classification_sig5`,
         disturbed0 = `1_classification_const`,
         disturbed1 = `1_classification_sig1`,
         disturbed2 = `1_classification_sig2`,
         disturbed3 = `1_classification_sig3`,
         disturbed4 = `1_classification_sig4`,
         disturbed5 = `1_classification_sig5`,
         bare0 = `2_classification_const` + `3_classification_const` + `4_classification_const`,
         bare1 = `2_classification_sig1` + `3_classification_sig1` + `4_classification_sig1`,
         bare2 = `2_classification_sig2` + `3_classification_sig2` + `4_classification_sig2`,
         bare3 = `2_classification_sig3` + `3_classification_sig3` + `4_classification_sig3`,
         bare4 = `2_classification_sig4` + `3_classification_sig4` + `4_classification_sig4`,
         bare5 = `2_classification_sig5` + `3_classification_sig5` + `4_classification_sig5`,
         natural0 = `5_classification_const`,
         natural1 = `5_classification_sig1`,
         natural2 = `5_classification_sig2`,
         natural3 = `5_classification_sig3`,
         natural4 = `5_classification_sig4`,
         natural5 = `5_classification_sig5`,
         wet0 = `6_classification_const` + `7_classification_const`,
         wet1 = `6_classification_sig1` + `7_classification_sig1`,
         wet2 = `6_classification_sig2` + `7_classification_sig2`,
         wet3 = `6_classification_sig3` + `7_classification_sig3`,
         wet4 = `6_classification_sig4` + `7_classification_sig4`,
         wet5 = `6_classification_sig5` + `7_classification_sig5`,)



### Build models

# make list of ee data types for making modlists
# classes <- c('0_classification',
#              '1_classification',
#              '2_classification',
#              '3_classification',
#              '4_classification',
#              '5_classification',
#              '6_classification',
#              '7_classification')
#
# weightings <- c('_const',
#                 '_sig1',
#                 '_sig2',
#                 '_sig3',
#                 '_sig4',
#                 '_sig5')
#
# vars <- outer(weightings, classes, rpaste0) %>%
#   as_tibble(.name_repair = 'universal') %>%
#   unite(modvars, everything(), sep = ' + ')
options(max.print = 9999999)
### Spring ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data_grouped %>%
  filter(Season == 'Spring') %>%
  mutate_at(72:101, scale) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')

#### Arachnida ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Arachnida + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod1 <- lmer(log(Arachnida + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod2 <- lmer(log(Arachnida + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod3 <- lmer(log(Arachnida + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod4 <- lmer(log(Arachnida + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod5 <- lmer(log(Arachnida + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE)
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
dev.off()
plot(allEffects(best.mod, residuals = TRUE),
     main = 'Arachnida, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Arachnida + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa5',
              'bare5', 'disturbed5', 'natural5', 'wet5')
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
fmod.red <- lmer(form.red, data = fD_spring_sub, REML = 'FALSE')
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
summary(ms1)
ms1
# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Arachnida + 1) ~ natural5 + (1|Site),
                data = fD_spring_sub, REML = FALSE)
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
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod1 <- lmer(log(Coccinellidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod2 <- lmer(log(Coccinellidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod3 <- lmer(log(Coccinellidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod4 <- lmer(log(Coccinellidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE,
             na.action = 'na.fail')
mod5 <- lmer(log(Coccinellidae + 1) ~
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
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)
png('spring_cocc_effect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(allEffects(best.mod, residuals = TRUE),
     main = '', #Coccinellidae - best landcover model
     partial.residual = list(lwd = 0),
     axes = list(x = list(natural4 =
                            list(lab =
                                   list(label =
                                          'Weighted proportion of \"natural\" landcover',
                                        cex = 1.5))),
                 y = list(lab = list(label = 'log(Coccinellidae density)', cex = 1.5))))
dev.off()
# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

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
ggsave('spring_cocc_groupvarweights.jpg', width = 7, height = 5)
# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Coccinellidae + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa4',
              'bare4', 'disturbed4', 'natural4', 'wet4')
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
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

# Build final model based on variable selection exercise.
# Best model is null. 2nd best is this one.
fin.mod <- lmer(log(Coccinellidae + 1) ~ natural4 + total_cover + (1|Site),
                                   data = fD_spring_sub, REML = FALSE)
summary(fin.mod)
png('spring_cocc_natEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('natural4', fin.mod, residuals = TRUE),
     main = 'Coccinellidae - best landcover + margin vegetation model',
     xlab = 'Weighted proportion of \"natural\" landcover',
     ylab = 'log(Coccinellidae density)')
dev.off()

png('spring_cocc_coverEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('total_cover', fin.mod, residuals = TRUE),
     main = '', #Coccinellidae - best landcover model
     partial.residual = list(lwd = 0),
     axes = list(x = list(total_cover =
                            list(lab =
                                   list(label =
                                          'Mean % plant cover in field margins',
                                        cex = 1.5))),
                 y = list(lab = list(label = 'log(Coccinellidae density)', cex = 1.5))))
dev.off()

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Ichneumonidae ####
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod1 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod2 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod3 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod4 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod5 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE)
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
#      main = 'Ichneumonidae, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Create global model that includes margin data.
# Global model is rank-deficient. We will have to 'trick' dredge.

# mod <- lmer(log(Ichneumonidae + 1) ~ shan + rich + total_cover +
#                     alfalfa4 + bare4 + disturbed4 + natural4 + wet4 +
#                     (1|Site),
#                  data = fD_spring_sub) # This doesn't work.

# List variables to include in global model.
vars.all <- c('shan', 'rich', 'total_cover', 'alfalfa0',
              'bare0', 'disturbed0', 'natural0', 'wet0')
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
fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 4), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Ichneumonidae + 1) ~
                  shan + (1|Site),
                data = fD_spring_sub)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Ichneumonidae, final model')

# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

#### Geocoris ####
options(na.action = 'na.fail')
# Build global models, keeping distanceWeights separate.
mod0 <- lmer(log(Geocoris + 1) ~
               alfalfa0 + bare0 + disturbed0 + natural0 + wet0 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod1 <- lmer(log(Geocoris + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod2 <- lmer(log(Geocoris + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod3 <- lmer(log(Geocoris + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod4 <- lmer(log(Geocoris + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod5 <- lmer(log(Geocoris + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

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
fmod.red <- lmer(form.red, data = fD_spring_sub, REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')
options(na.action = 'na.omit')
# Build final model based on variable selection exercise.
fin.mod <- lmer(log(Geocoris + 1) ~ natural2 + total_cover +
                  (1|Site),
                data = fD_spring, REML = FALSE)
summary(fin.mod)
plot(allEffects(fin.mod, residuals = TRUE), main = 'Geocoris, final model')

png('spring_geo_natEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('natural2', fin.mod, residuals = TRUE),
     main = NULL,
     xlab = 'Weighted proportion of \"natural\" landcover',
     ylab = 'log(Geocoris density)')
dev.off()
png('spring_geo_coverEffect.jpg', width = 7, height = 5, units = 'in', res = 300)
plot(effect('total_cover', fin.mod, residuals = TRUE),
     main = '', #Coccinellidae - best landcover model
     partial.residual = list(lwd = 0),
     axes = list(x = list(total_cover =
                            list(lab =
                                   list(label =
                                          'Mean % cover in field margins',
                                        cex = 1.5))),
                 y = list(lab = list(label = 'log(Geocoris density)', cex = 1.5))))
dev.off()
# Clean environment before moving on to next taxon.
rm(mod0, mod1, mod2, mod3, mod4, mod5, cand_mods, dredges, mod_table, best.mod,
   importance_tab, p, group_importance, vars.all, form, vars.red, form.red,
   fmod.red, ms1, fin.mod)

### Fall ####

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
             REML = FALSE)
mod1 <- lmer(log(Arachnida + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod2 <- lmer(log(Arachnida + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod3 <- lmer(log(Arachnida + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod4 <- lmer(log(Arachnida + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod5 <- lmer(log(Arachnida + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

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
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE)
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
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
             REML = FALSE)
mod1 <- lmer(log(Coccinellidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod2 <- lmer(log(Coccinellidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod3 <- lmer(log(Coccinellidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod4 <- lmer(log(Coccinellidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod5 <- lmer(log(Coccinellidae + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

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
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE)
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
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
             REML = FALSE)
mod1 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod2 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod3 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod4 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod5 <- lmer(log(Ichneumonidae + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

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
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
             data = fD_fall,               REML = FALSE)
mod1 <- lmer(log(Geocoris + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,               REML = FALSE)
mod2 <- lmer(log(Geocoris + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,               REML = FALSE)
mod3 <- lmer(log(Geocoris + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,               REML = FALSE)
mod4 <- lmer(log(Geocoris + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,               REML = FALSE)
mod5 <- lmer(log(Geocoris + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,               REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

p <- ggplot(data = group_importance,
            aes(x = '', y = distWeight, fill = weight)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

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
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
             REML = FALSE)
mod1 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod2 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod3 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod4 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_spring,
             REML = FALSE)
mod5 <- lmer(log(Acyrthosiphon + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_spring,
             REML = FALSE)
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
mod_table

# Extract and examine the best model.
best.mod <- get.models(mod_table, subset = 1)[[1]]
summary(best.mod)

# plot(allEffects(best.mod, residuals = TRUE),
#      main = 'Acyrthosiphon, best landcover model')

# Plot the variable importance.
importance_tab <- sw(mod_table) %>%
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

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
                 REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 5), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
             REML = FALSE)
mod1 <- lmer(log(NonAcy + 1) ~
               alfalfa1 + bare1 + disturbed1 + natural1 + wet1 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod2 <- lmer(log(NonAcy + 1) ~
               alfalfa2 + bare2 + disturbed2 + natural2 + wet2 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod3 <- lmer(log(NonAcy + 1) ~
               alfalfa3 + bare3 + disturbed3 + natural3 + wet3 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod4 <- lmer(log(NonAcy + 1) ~
               alfalfa4 + bare4 + disturbed4 + natural4 + wet4 + (1|Site),
             data = fD_fall,
             REML = FALSE)
mod5 <- lmer(log(NonAcy + 1) ~
               alfalfa5 + bare5 + disturbed5 + natural5 + wet5 + (1|Site),
             data = fD_fall,
             REML = FALSE)
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
  tidy() %>%
  arrange(names) %>%
  separate(names, c('class', 'distWeight'), sep = -1) %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `0` = 'constant',
                                       `1` = 'aggressive',
                                       `2` = 'moderately aggressive',
                                       `3` = 'moderate',
                                       `4` = 'slight',
                                       `5` = 'minimal')))

p <- ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

# Plot the variable importance, summarized by distanceWeight.
group_importance <- importance_tab %>%
  group_by(distWeight) %>%
  summarize(weight = sum(x))

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
fmod.red <- lmer(form.red, data = fD_fall_sub, REML = FALSE)
# Replace reduced model formula with full global model formula.
attr(fmod.red@frame, "formula") <- form
# Check formula attribute of reduced model.
formula(fmod.red) # Looks good.

# Run dredge() with m.max parameter to avoid convergence failures.
ms1 <- model.sel(lapply(dredge(fmod.red, evaluate = FALSE, m.max = 2), eval))
summary(ms1)

# Plot the variable importance.
importance_tab <- sw(ms1) %>%
  tidy() %>%
  mutate(names = str_replace(names, '[:digit:]', '')) %>%
  rename(varWeight = x)

p <- ggplot(data = importance_tab,
            aes(x = reorder(names,-varWeight), y = varWeight)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) +
  labs(x = 'Variable', y = 'Importance')

ggplotly(p, tooltip = 'x')

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
  tidy()

p <- ggplot(data = importance_tab, aes(x = reorder(names, -x), y = '', fill = x)) +
  geom_tile() +
  theme(axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(angle = 45)) +
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       fill = 'Variable importance')

ggplotly(p, tooltip = 'x')

best.delta.mod = lmer(AllAph ~
                         Ichneumonidae + (1|Site),
                      data = delta_density,
                      na.action = 'na.fail',
                      REML = FALSE)

plot(allEffects(best.delta.mod, residuals = TRUE))
summary(best.delta.mod)
