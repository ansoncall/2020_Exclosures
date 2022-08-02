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
# landcover <- read_csv('tidy_data/landcover.csv', col_types = 'ffffdddddddddddd')
landcover <- read_csv('tidy_data/landcover_fixed.csv', #######################
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
             'weedy',
             'wet',
             'water')


# Collapse classes ########################################################
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

### fall data only ####

fD_fall <- field_data %>%
  filter(Season == 'Fall') %>%
  mutate_at(21:69, ~ as.vector(scale(.))) # All landcover scores are scaled here.

# Subset data for modeling with landcover classes and margin data.
fD_fall_sub <- fD_fall %>% filter(Site != 'Yerington')

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

# build mod tabs

# Add models with margin data to model selection table.
makeGlobalLandcoverModTab <- function(mod_table, rank = 1){
  # mod_table = arachnida_mod_table_sp
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






# list of taxa to build tables for
taxa <- c('Acyrthosiphon', 'AllAph', 'Arachnida', 'Coccinellidae', 'Geocoris')

# list to fill with modTables
noMarginTabs <- list()

# build spring tabs
for (i in 1:length(taxa)) {

  noMarginTabs[[i]] <- buildLandcoverModTab(taxa[[i]], fD_spring, 4)

}

# build fall tabs
for (i in 1:length(taxa)) {

  noMarginTabs[[i]] <- buildLandcoverModTab(taxa[[i]], fD_fall, 4)

}

# currently: 8 class, fixed
itemNames <- paste(taxa, '8class_fixed', sep = "_") ################################
for (i in 1:length(taxa)) {

  saveRDS(noMarginTabs[[i]], itemNames[[i]])

}
