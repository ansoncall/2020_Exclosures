# This script is for generating and saving model selection tables for landcover
# and field margin data. This makes the data_analysis script easier to tinker
# with, because you don't need to refit the models every time you make a change.

library(car) # for Anova() on lmer model objects
library(crayon) # for colored terminal outputs
library(hardhat) # for get_levels to extract factor levels in a tidy way
library(lme4) # for univariate mixed-effects models
library(lmerTest) # for lmer with p values
library(magrittr) # for assignment and exposition pipes
library(MuMIn) # model selection tools
library(mvabund) # for building multivariate mods of insect density
library(plotly) # interactive plots with plotly()
library(sjPlot) # create effects plots on lmer objects
library(tidyverse) # R packages for data science


# Note - I'm doing this the dumb way (changing vars by hand instead of looping)
# because I don't have time to write new code
# look for ##############################################33333

# read in data ####
data <- read_csv('tidy_data/data.csv', col_types = 'fffffdddddddddddddddddddd')
data_long <- read_csv('tidy_data/data_long.csv', col_types = 'fffffffd')
landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffdddddddddddd')
landcoverFixed <- read_csv('tidy_data/landcover_fixed.csv', #######################
                      col_types = 'ffffdddddddddddd')
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

# define varlist
varList <- c('alfalfa', ##################################################
             'naturalArid',
             'dirt',
             'ag',
             'impermeable',
             'weedyWet',
             'water')


# Collapse classes ########################################################
landcover %<>%
  # 2022-08-03T20:15:27Z Fix error: distWeights 1 and 2 mixed up
  mutate(distanceWeight = case_when(distanceWeight == 'sig1' ~ 'sig2',
                                    distanceWeight == 'sig2' ~ 'sig1',
                                    distanceWeight == 'no' ~ 'no',
                                    distanceWeight == 'const' ~ 'const',
                                    distanceWeight == 'sig3' ~ 'sig3',
                                    distanceWeight == 'sig4' ~ 'sig4',
                                    distanceWeight == 'sig5' ~ 'sig5')) %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedyWet = class5 + class8,
         water = class11
  ) %>%
  select(-(class0:class11)) %>%
  mutate(distanceWeight = fct_relevel(distanceWeight, 'no', after = Inf))
levels(landcover$distanceWeight)

landcoverFixed %<>%
  # 2022-08-03T20:15:27Z Fix error: distWeights 1 and 2 mixed up
  mutate(distanceWeight = case_when(distanceWeight == 'sig1' ~ 'sig2',
                                    distanceWeight == 'sig2' ~ 'sig1',
                                    distanceWeight == 'no' ~ 'no',
                                    distanceWeight == 'const' ~ 'const',
                                    distanceWeight == 'sig3' ~ 'sig3',
                                    distanceWeight == 'sig4' ~ 'sig4',
                                    distanceWeight == 'sig5' ~ 'sig5')) %>%
  mutate(alfalfa = class6,
         naturalArid = class10,
         dirt = class2 + class2 + class7,
         ag = class0 + class3,
         impermeable = class4 + class9,
         weedyWet = class5 + class8,
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

landcover_longFixed <- landcoverFixed %>%
  pivot_longer(alfalfa:water,
               names_to = 'class',
               values_to = 'areaScore') %>%
  select(siteId, distanceWeight, site, field, class, areaScore)


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

landcover_wideFixed <- landcoverFixed %>%
  pivot_wider(id_cols = c(site, field),
              names_from = distanceWeight,
              values_from = alfalfa:water) %>%
  mutate(id = paste0(site, field), .before = site)
# join landcover data
field_data <- left_join(mean_density_field, landcover_wide)
field_dataFixed <- left_join(mean_density_field, landcover_wideFixed)


### spring data only ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

fD_springFixed <- field_dataFixed %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')
fD_spring_subFixed <- fD_springFixed %>% filter(Site != 'Yerington')

### fall data only ####

fD_fall <- field_data %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

fD_fallFixed <- field_dataFixed %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')
fD_fall_subFixed <- fD_fallFixed %>% filter(Site != 'Yerington')

# functions ####
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



# list of taxa to build tables for
taxa <- c('Acyrthosiphon', 'NonAcy', 'AllAph', 'Arachnida', 'Coccinellidae',
          'Geocoris', 'Ichneumonidae')



# spring List
springNonFix <- list()
springItemNames <- paste(taxa, 'spring_7class', sep = "_")
springFix <- list()
springFixItemNames <- paste(taxa, 'spring_7class_fixed', sep = "_")
# build spring non fixed tabs
for (i in 1:length(taxa)) {

  springNonFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_spring, 3)

}


# build spring fixed tabs
for (i in 1:length(taxa)) {

  springFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_springFixed, 3)
}

# fall List
fallNonFix <- list()
fallItemNames <- paste(taxa, 'fall_7class', sep = "_")
fallFix <- list()
fallFixItemNames <- paste(taxa, 'fall_7class_fixed', sep = "_")
# build fall non fixed tabs
for (i in 1:length(taxa)) {

  fallNonFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_fall, 3)

}


# build fall fixed tabs
for (i in 1:length(taxa)) {

  fallFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_fallFixed, 3)

}

tabsList <- c(springNonFix, springFix, fallNonFix, fallFix)
namesList <- c(springItemNames, springFixItemNames, fallItemNames, fallFixItemNames)
# currently: 7 class
 ################################


for (i in 1:length(namesList)) {

  saveRDS(tabsList[[i]], paste0('modTabs-distFixed/',namesList[[i]]))

}



  # 8 class here ####
# read in data ####
data <- read_csv('tidy_data/data.csv', col_types = 'fffffdddddddddddddddddddd')
data_long <- read_csv('tidy_data/data_long.csv', col_types = 'fffffffd')
landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffdddddddddddd')
landcoverFixed <- read_csv('tidy_data/landcover_fixed.csv', #######################
                           col_types = 'ffffdddddddddddd')
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

mean_density <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Plot, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  separate(id, c('Site', 'Field', 'Plot', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density)

# define varlist
varList <- c('alfalfa', ##################################################
             'naturalArid',
             'dirt',
             'ag',
             'impermeable',
             'weedy',
             'wet',
             'water')


# Collapse classes ########################################################
landcover %<>%
  # 2022-08-03T20:15:27Z Fix error: distWeights 1 and 2 mixed up
  mutate(distanceWeight = case_when(distanceWeight == 'sig1' ~ 'sig2',
                                    distanceWeight == 'sig2' ~ 'sig1',
                                    distanceWeight == 'no' ~ 'no',
                                    distanceWeight == 'const' ~ 'const',
                                    distanceWeight == 'sig3' ~ 'sig3',
                                    distanceWeight == 'sig4' ~ 'sig4',
                                    distanceWeight == 'sig5' ~ 'sig5')) %>%
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

landcoverFixed %<>%
  # 2022-08-03T20:15:27Z Fix error: distWeights 1 and 2 mixed up
  mutate(distanceWeight = case_when(distanceWeight == 'sig1' ~ 'sig2',
                                    distanceWeight == 'sig2' ~ 'sig1',
                                    distanceWeight == 'no' ~ 'no',
                                    distanceWeight == 'const' ~ 'const',
                                    distanceWeight == 'sig3' ~ 'sig3',
                                    distanceWeight == 'sig4' ~ 'sig4',
                                    distanceWeight == 'sig5' ~ 'sig5')) %>%
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

landcover_longFixed <- landcoverFixed %>%
  pivot_longer(alfalfa:water,
               names_to = 'class',
               values_to = 'areaScore') %>%
  select(siteId, distanceWeight, site, field, class, areaScore)


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

landcover_wideFixed <- landcoverFixed %>%
  pivot_wider(id_cols = c(site, field),
              names_from = distanceWeight,
              values_from = alfalfa:water) %>%
  mutate(id = paste0(site, field), .before = site)
# join landcover data
field_data <- left_join(mean_density_field, landcover_wide)
field_dataFixed <- left_join(mean_density_field, landcover_wideFixed)


### spring data only ####

# Subset data for modeling with landcover classes.
fD_spring <- field_data %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

fD_springFixed <- field_dataFixed %>%
  filter(Season == 'Spring') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_spring_sub <- fD_spring %>% filter(Site != 'Yerington')
fD_spring_subFixed <- fD_springFixed %>% filter(Site != 'Yerington')

### fall data only ####

fD_fall <- field_data %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

fD_fallFixed <- field_dataFixed %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.
#2022-08-10T15:48:14Z
## ERROR ####
# surface water class not scaled here??? not fixed yet!!

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')
fD_fall_subFixed <- fD_fallFixed %>% filter(Site != 'Yerington')

# functions ####
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


# list of taxa to build tables for
taxa <- c('Acyrthosiphon', 'NonAcy', 'AllAph', 'Arachnida', 'Coccinellidae',
          'Geocoris', 'Ichneumonidae')



# spring List
springNonFix <- list()
springItemNames <- paste(taxa, 'spring_8class', sep = "_")
springFix <- list()
springFixItemNames <- paste(taxa, 'spring_8class_fixed', sep = "_")
# build spring non fixed tabs
for (i in 1:length(taxa)) {

  springNonFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_spring, 3)

}


# # build spring fixed tabs
# for (i in 1:length(taxa)) {
#
#   springFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_springFixed, 3)
#
# }

# fall List
fallNonFix <- list()
fallItemNames <- paste(taxa, 'fall_8class', sep = "_")
fallFix <- list()
fallFixItemNames <- paste(taxa, 'fall_8class_fixed', sep = "_")
# build fall non fixed tabs
for (i in 1:length(taxa)) {

  fallNonFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_fall, 3)

}


# # build fall fixed tabs
# for (i in 1:length(taxa)) {
#
#   fallFix[[i]] <- buildLandcoverModTab(taxa[[i]], fD_fallFixed, 3)
#
# }
#
# tabsList <- c(springNonFix, springFix, fallNonFix, fallFix)
tabsList <- c(springNonFix, fallNonFix)
# namesList <- c(springItemNames, springFixItemNames, fallItemNames, fallFixItemNames)
namesList <- c(springItemNames, fallItemNames)
# # currently: 8 class
# ################################
#
#
for (i in 1:length(namesList)) {

  saveRDS(tabsList[[i]], paste0('modTabs/',namesList[[i]]))

}
