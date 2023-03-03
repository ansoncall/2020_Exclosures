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
landcoverFixed <- read_csv('tidy_data/landcover_fixed.csv',
                      col_types = 'ffffdddddddddddd')
veg_plots <- read_csv('tidy_data/veg_plots.csv', col_types = 'fffffff')
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

# prepare landcover data ####
mean_density_field <- data_long %>%
  mutate(Density = case_when(Treatment =='Pre-' ~ Density/3,
                             Treatment !='Pre-' ~ Density)) %>%
  unite(id, Site, Field, Season, Taxa, remove = FALSE) %>%
  group_by(id) %>%
  summarize(Mean_Density = mean(Density)) %>%
  separate(id, c('Site', 'Field', 'Season', 'Taxa')) %>%
  pivot_wider(names_from = Taxa, values_from = Mean_Density) %>%
  mutate(id = paste0(Site, 0, Field), .before = Site)

# make list of input landcover data tables
lcList <- list(landcover, landcoverFixed)
names(lcList) <- c('landcover', 'landcoverFixed')
# list to hold outputs
lcListOut <- list()
# must include 7-class and 8-class version for each input
# 7 class = weedy+wet combined; 8 class = weedy+wet separate
# define varlists
# 7 class
varList7 <- c('alfalfa', ##################################################
              'naturalArid',
              'dirt',
              'ag',
              'impermeable',
              'weedyWet',
              'water')
# 8 class
varList8 <- c('alfalfa', ##################################################
              'naturalArid',
              'dirt',
              'ag',
              'impermeable',
              'weedy',
              'wet',
              'water')
# make list of varLists
varLists <- list(varList7, varList8)
names(varLists) <- c('varList7', 'varList8')

# loop over landcover inputs
for (i in 1:length(lcList)) {
  print(names(lcList[i]))

  landcover <- lcList[[i]]

  landcover7 <- landcover %>%
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

  landcover8 <- landcover %>%
    mutate(alfalfa = class6,
           naturalArid = class10,
           dirt = class2 + class2 + class7,
           ag = class0 + class3,
           impermeable = class4 + class9,
           weedy =  class5,
           wet = class8,
           water = class11
    ) %>%
    select(-(class0:class11)) %>%
    mutate(distanceWeight = fct_relevel(distanceWeight, 'no', after = Inf))
  levels(landcover$distanceWeight)

  # collect 7-class and 8-class data into a list
  lcClassList <- list(landcover7, landcover8)
  lcClassListOut <- list()
  # loop over 7-class and 8-class data
  for (j in 1:length(lcClassList)) {

    landcover <- lcClassList[[j]]
    # reshape
    # lengthen
    landcover_long <- landcover %>%
      pivot_longer(alfalfa:water,
                   names_to = 'class',
                   values_to = 'areaScore') %>%
      select(siteId, distanceWeight, site, field, class, areaScore)
    # widen
    landcover_wide <- landcover %>%
      pivot_wider(id_cols = c(site, field),
                  names_from = distanceWeight,
                  values_from = alfalfa:water) %>%
      mutate(id = paste0(site, field), .before = site)

    # join landcover data to arthropod data
    field_data <- left_join(mean_density_field, landcover_wide)

    lcClassListOut[[j]] <- field_data
  }
  lcListOut[[i]] <- lcClassListOut
}

# flatten lists and add names
landCoverList <- unlist(lcListOut, recursive=FALSE)
names(landCoverList) <- c('regular7', 'regular8', 'fixed7', 'fixed8')

rm(lcClassList, lcClassListOut, lcList, lcListOut)

# scale all landcover vars in all inputs
landCoverScaled <- list()
for (i in 1:length(landCoverList)) {
  landCoverScaled[[i]] <- landCoverList[[i]] %>%
    mutate(across(contains('_'), ~ as.vector(scale(.))))
}

# create '_sub' (no Yerington) datasets
landCoverScaledSub <- list()
for (i in 1:length(landCoverScaled)) {
  landCoverScaledSub[[i]] <-
    landCoverScaled[[i]] %>% filter(Site != 'Yerington')
}

# join 'Scaled' and 'ScaledSub'
allScaled <- c(landCoverScaled, landCoverScaledSub)
# name all elements
names(allScaled) <- paste0(rep(names(landCoverList), 2),
                        c(rep('_full', 4), rep('_sub', 4)))
# # check
# View(allScaled)

# split into spring and fall data
allScaledSpring <- list()
for (i in 1:length(allScaled)) {
  allScaledSpring[[i]] <- allScaled[[i]] %>%
    filter(Season == 'Spring')
}

allScaledFall <- list()
for (i in 1:length(allScaled)) {
  allScaledFall[[i]] <- allScaled[[i]] %>%
    filter(Season == 'Fall')
}

allLandCover <- c(allScaledSpring, allScaledFall)
names(allLandCover) <- paste0(rep(names(allScaled), 2),
                              c(rep('_spring', 8), rep('_fall', 8)))
# # check
# View(allLandCover)

# remove helper lists
rm(allScaled, allScaledFall, allScaledSpring,
   landCoverList, landCoverScaled, landCoverScaledSub,
   landcover7, landcover8, landcover_long, landcover_wide,
   landcover, mean_density, mean_density_field)
i = 1
rm(taxon,data,m.max,i)
# functions ####
buildLandcoverModTab <- function(taxon = 'empty', data = 'empty', m.max = 3){
  # taxon='NonAcy'
  # data=allLandCover[[i]]
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
          'Geocoris', 'Ichneumonoidea')


modLoop <- list()
for (i in 1:length(allLandCover)) {

  innerModLoop <- list()
  data <- allLandCover[[i]]
  for (j in 1:length(taxa)) {

    innerModLoop[[j]] <- buildLandcoverModTab(taxa[[j]], data, 3)

  }

  modLoop[[i]] <- innerModLoop
}

# get the names right!
names(modLoop) <- names(allLandCover)

for (i in 1:length(modLoop)) {

  names(modLoop[[i]]) <- taxa

}

View(modLoop)

# export RDS
for (i in 1:length(modLoop)) {

  dataset <- names(modLoop[i])

  for (j in 1:length(modLoop[[i]])) {

    modTabList <- modLoop[[i]]
    response <- names(modTabList[j])

    saveRDS(modTabList[[j]], paste0('modTabs/',
                                    response,
                                    '_',
                                    dataset))

  }

}
