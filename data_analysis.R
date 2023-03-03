# this is the main analysis script. it depends on the data_processing script for
# tidy data inputs, and the data_analysis_landcover_tables script for
# pre-generated AICc tables of arthropod~landcover models.

# Load packages ####
library(BAMMtools) # for fast Jenks breaks
library(broom.mixed) # for tidying model outputs
library(car) # for Anova() on lmer model objects
library(classInt) # for kmeans clustering
library(corrr) # for correlation plots of landcover vars
library(crayon) # for colored terminal outputs
library(data.table) # for rbindlist() to rbind a list of tables and make id col
library(DiagrammeR) # for SEM plots
library(DHARMa) # for simulated residuals
library(dotwhisker) # for effect size plots
library(effects) # for effects plots
library(emmeans) # for computing SEM marginal means
library(ggeffects) # for easy effects plots
library(ggfortify) # create PCA plots
library(ggiraphExtra) # more easy effects plots
library(glmmTMB) # for mixed-effects GLM
library(ggpmisc) # make ggplots that show R2 value with stat_poly_eq()
library(ggridges) # for ggridges plots
library(grid) # for grobTree() to make text annotations on violin plots
library(gridExtra) # create multi-panel plots
library(gtools) # for mixedsort() to arrange factor levels in vegdata tibble
library(hardhat) # for get_levels to extract factor levels in a tidy way
library(knitr) # for knitting R markdown docs
library(lavaan) # for path analysis/SEM
library(lavaanPlot) # for plotting lavaan models
library(lme4) # for univariate mixed-effects models
library(lmerTest) # for lmer with p values
library(magrittr) # for assignment and exposition pipes
library(MuMIn) # model selection tools
library(mvabund) # for building multivariate mods of insect density
library(parallel) # parallelization of dredge()
library(performance) # for pseudo-R^2 of GLMERs
library(plotly) # interactive plots with plotly()
library(sjmisc) # for to_dummy function
library(sjPlot) # create effects plots on lmer objects, plot tables
library(tidyselect) # for peek()
library(tidytable) # for get_dummies to make SEM data
library(tidytext) # sort ggcols after faceting
library(tidyverse) # R packages for data science
library(varhandle) # easily create dummy vars with to.dummy()
library(vegan) # for diversity indices in vegdata
library(webshot) # to capture tab_model output (or other html output) as png


# Define functions ####
# note: this is not the only place functions are defined

# p value formatting function
# not sure how this works. from stackexchange
pvalr <- function(pvals, sig_limit = .001, digits = 3, html = FALSE) {

  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0("%.", digits, "f"), x)
    zzz <- paste0("0.", paste(rep("0", digits), collapse = ""))
    res[res == paste0("-", zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig_limit) {
    if (x < sig_limit)
      if (html)
        return(sprintf("&lt; %s", format(sig_limit))) else
          return(sprintf("< %s", format(sig_limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig_limit = sig_limit)
}

# define stupid "reverse paste" function for mapping paste0() over a list
rpaste0 <- function(x, y) {
  paste0(y, x)
}

# import/wrangle data ####
# subplot-level data, all vars
# pre- data not /3, can use area offset to correct.
# No transformations have been applied
subplot_data_raw <- read_csv("tidy_data/subplotDataRaw.csv",
                           col_types = "ffffff")

# "Pre-" density has already been /3
# No transformations have been applied
subplot_data <- read_csv("tidy_data/subplotData.csv",
                        col_types = "ffffff")

# make plot-level data, excluding "Pre-" measurements
plot_data <- subplot_data %>%
  filter(Treatment != "Pre-") %>%
  group_by(Site, Field, Plot, Season) %>%
  summarize(across(where(is.numeric), .fns = mean), .groups = "keep")

# make field-level data by taking means across fields (all plots pooled)
field_data <- subplot_data %>%
  filter(Treatment != "Pre-") %>%
  group_by(Site, Field, Treatment, Season) %>%
  summarize(across(where(is.numeric), .fns = mean), .groups = "keep")

# raw data from vegetation plots
veg_plots <- read_csv("tidy_data/vegPlots.csv")


# define color palette ####
# Acyrthosiphon aphids: #548235
# Non-Acyrthosiphon aphids: #3B3838
# Spring: #76db91
# Fall: #9e3c21
# Predators: "Spectral"
# Sites: "Set1"

# Data exploration ####
## Arthropod data summary
source("data_analysis_arthropod_exploration.R", echo = TRUE)

## Landcover data summary
source("data_analysis_classification_summary.R", echo = TRUE)



# Wrangle data for model selection ####
# split spring and fall data
df_sp <- subplot_data_raw %>%
  # spring only
  filter(Season == "Spring") %>%
  # "regular" landcover only
  select(!contains("fix")) %>%
  # log AllAph col
  mutate(log_AllAph = log(AllAph + 1)) %>%
  # center and scale (not needed with rank transform)
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

df_fa <- subplot_data_raw %>%

  filter(Season == "Fall") %>%
  # "regular" landcover only
  select(!contains("fix")) %>%
  # log AllAph col
  mutate(log_AllAph = log(AllAph + 1)) %>%
  # center and scale (not needed with rank transform)
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

# use rank transform on landcover vars
df_sp_rnk <- subplot_data_raw %>%
  # spring only
  filter(Season == "Spring") %>%
  # "regular" landcover only
  select(!contains("fix")) %>%
  # rank transform landcover to uniformly distribute
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~dense_rank(.x))) %>%
  # log AllAph col
  mutate(log_AllAph = log(AllAph + 1)) %>%
  # # center and scale (not needed with rank transform)
  # mutate(across(.cols = contains("_"), # all landcover + log_AllAph
  #               .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

df_fa_rnk <- subplot_data_raw %>%
  # spring only
  filter(Season == "Fall") %>%
  # "regular" landcover only
  select(!contains("fix")) %>%
  # rank transform landcover to uniformly distribute
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~dense_rank(.x))) %>%
  # log AllAph col
  mutate(log_AllAph = log(AllAph + 1)) %>%
  # # center and scale (not needed with rank transform)
  # mutate(across(.cols = contains("_"), # all landcover + log_AllAph
  #               .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

# # example dotcharts - shows distribution of explanatory variables
# dotchart(sort(df_sp_rnk$AllAph))
# dotchart(sort(df_sp$log_AllAph))
# dotchart(sort(df_fa$log_AllAph))

# Fit models ####

rebuild <- askYesNo("Would you like to rebuild model selection tables?")

if (rebuild == TRUE) {

  ## source external scripts to build model selection tables
  ## set number of cores to be used in parallel processing
  n_cores <- detectCores() - 4
  # glmer, negative binomial, scaled vars
  source("nbMixed.R", echo = TRUE) # system.time 161 seconds
  # glmer, negative binomial, ranked vars
  source("nbMixedRanked.R", echo = TRUE)
  # glmer, poisson, scaled vars
  source("poisMixed.R", echo = TRUE)
  # glmer, poisson, ranked vars
  source("poisMixedRanked.R", echo = TRUE)

  ## SOURCE ####
  # nb mixed mods, landcover scaled
  source("aphNbMixed.R", echo = TRUE)

  save.image(file = "analysis_env.RData")

} else {

  # Optional: start here ####
  load("analysis_env.RData")
}

# Collect top predator models ####

# collects the top models from each model family and compares them in a new
# model selection table.
source("collectMods_preds.R", echo = TRUE)

# # Compare top predator models
# anth_fams_sp %>% View # nb_scaled by at least delta>2
# anth_fams_fa %>% View # nb_scaled by delta 1.27. NO RANDOM EFFECT in top mod
# ara_fams_sp %>% View # both pois mods close, and they disagree
# ara_fams_fa %>% View # nb mods agree, pois mods are delta+ 15
# cocc_fams_sp %>% View
        # scaled mods agree, ranked mods differ but delta +4 anyway
# cocc_fams_fa %>% View # scaled mods agree, ranked mods differ,
#                      # ranked have slightly better fit but deltas are close
# geo_fams_sp %>% View # nb_scaled by delta+ 13
# geo_fams_fa %>% View # all mods agree and are generally close
# ich_fams_sp %>% View # mods mostly agree and deltas are close
# ich_fams_fa %>% View # nb_scaled by delta+7 NO RANDOM EFFECT in top mod

# Make table of top (no veg) predator models ####
# make list of best models
best_mod_list <- list(
  "best.ant.sp" = get.models(nb_scaled$tab_nb_anth_sp_scaled, 1)[[1]],
  "best.ara.sp" = get.models(nb_scaled$tab_nb_ara_sp_scaled, 1)[[1]],
  "best.coc.sp" = get.models(nb_scaled$tab_nb_cocc_sp_scaled, 1)[[1]],
  "best.geo.sp" = get.models(nb_scaled$tab_nb_geo_sp_scaled, 1)[[1]],
  "best.ich.sp" = get.models(nb_scaled$tab_nb_ich_sp_scaled, 1)[[1]],
  "best.ant.fa" = get.models(nb_scaled$tab_nb_anth_fa_scaled, 1)[[1]],
  "best.ara.fa" = get.models(nb_scaled$tab_nb_ara_fa_scaled, 1)[[1]],
  "best.coc.fa" = get.models(nb_scaled$tab_nb_cocc_fa_scaled, 1)[[1]],
  "best.geo.fa" = get.models(nb_scaled$tab_nb_geo_fa_scaled, 1)[[1]],
  "best.ich.fa" = get.models(nb_scaled$tab_nb_ich_fa_scaled, 1)[[1]]
)

# build empty tibble to hold stats
stats_df <- tibble(Taxon = rep(c("Anthocoridae", "Arachnida", "Coccinellidae",
                                "Geocoris", "Ichneumodoidea"), 2),
                  Season = c(rep("Spring", 5), rep("Fall", 5)),
                  MarginalR2 = c(0),
                  ConditionalR2 = c(0),
                  effects1 = c("none"),
                  effects2 = c("none"),
                  coefs1 = c(0),
                  coefs2 = c(0))

# fill tibble with stats
for (i in seq_along(best_mod_list)){
  stats_df$MarginalR2[[i]] <- r2(best_mod_list[[i]])[[2]]
  stats_df$ConditionalR2[[i]] <- r2(best_mod_list[[i]])[[1]]
  stats_df$effects1[[i]] <- names(best_mod_list[[i]]$frame)[2]
  stats_df$effects2[[i]] <- names(best_mod_list[[i]]$frame)[3]
  stats_df$coefs1[[i]] <- fixef(best_mod_list[[i]])$cond[2]
  stats_df$coefs2[[i]] <- fixef(best_mod_list[[i]])$cond[3]
}


stats_df %>%
  group_by(Season) %>%
  tab_df(title = "Top predator models (no vegetation data included)")


# Review predator models ####

# best to review these by hand. change input models manually.

# # optional: review a single mod table
# nb_scaled$tab_nb_cocc_sp_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>% View

# choose model to review
review_mod <- get.models(nb_scaled$tab_nb_cocc_sp_scaled, 1)[[1]]

# show summary
summary(review_mod) # no random effect variance. essentially equivalent to
                         # a fixed effects mod. I checked.
# basic effects plots
plot(allEffects(review_mod, residuals = TRUE))

# tidy coeffs plot
tidy(review_mod, conf.int = TRUE) %>%
  # needs "model" column for dwplot
  mutate(model = 1) %>%
  dwplot


# extract residuals, fitted values
pearson_res <- resid(review_mod, type = "pearson")
working_res <- resid(review_mod, type = "working")
default_res <- resid(review_mod)
dharma_res <- simulateResiduals(review_mod)
fitted <- fitted(review_mod)

# basic DHARMa plot
plot(dharma_res)

# other residual plots
plot(fitted, pearson_res)
plot(fitted, default_res)


# Predator models with vegdata ####
# drop Yerington (NA vegdata values)
df_sp_vd <- df_sp %>% filter(!is.na(shan))
df_fa_vd <- df_fa %>% filter(!is.na(shan))
# # Check nrow
# df_sp %>% nrow
# df_sp_vd %>% nrow
# df_fa %>% nrow
# df_fa_vd %>% nrow

## Add vegdata to top mods
source("compareVeg_pred.R", echo = TRUE)

### Spring
# Anthocoridae
r2(best_mod_list$best.ant.sp)
r2(best.ant.sp.vd) ## new mod has one less factor
# vedict - keep original

# Arachnida
r2(best_mod_list$best.ara.sp)
r2(best.ara.sp.vd) ## same model structures. original has more data but worse r2
# vedict - keep original

# Coccinellidae
r2(best_mod_list$best.coc.sp)
r2(best.coc.sp.vd) ## same model structures. original has more data and better
                   ## marginal r2
# vedict - keep original

# Ichneumonoidea
r2(best_mod_list$best.ich.sp)
r2(best.ich.sp.vd) ## new model has one less factor. new has better marginal r2,
                   ## but worse conditional r2
# vedict - keep original

### Fall
# Anthocoridae
r2(best_mod_list$best.ant.fa)
r2(best.ant.fa.vd) ## total cover looking like the best predictor here
## verdict - OVERTURN. new model is better!

# Arachnida
# must drop wateringMethod here because all sites are flooded.
r2(best_mod_list$best.ara.fa)
r2(best.ara.fa.vd) ## new model with total cover instead of wateringMethod is
                   ## better.
## verdict - OVERTURN. new model is better!


# Coccinellidae
r2(best_mod_list$best.coc.fa)
r2(best.coc.fa.vd) ## new model has one less factor.
## verdict - keep original

# Ichneumonoidea
r2(best_mod_list$best.ich.fa)
r2(best.ich.fa.vd) ## same models. original has better marginal r2
## verdict - keep original

### Summarize best pred models (w/vegdata)
new_mods <- list("antFA" = best.ant.fa.vd, "araFA" = best.ara.fa.vd)

# build empty tibble to hold stats
stats_df_new <- tibble(Taxon = c("Anthocoridae", "Arachnida"),
                  Season = rep("Fall", 2),
                  MarginalR2 = c(0),
                  ConditionalR2 = c(0),
                  effects1 = c("none"),
                  effects2 = c("none"),
                  coefs1 = c(0),
                  coefs2 = c(0))

# fill tibble with stats
for (i in seq_along(new_mods)){
  stats_df_new$MarginalR2[[i]] <- r2(new_mods[[i]])[[2]]
  stats_df_new$ConditionalR2[[i]] <- r2(new_mods[[i]])[[1]]
  stats_df_new$effects1[[i]] <- names(new_mods[[i]]$frame)[2]
  stats_df_new$effects2[[i]] <- names(new_mods[[i]]$frame)[3]
  stats_df_new$coefs1[[i]] <- fixef(new_mods[[i]])$cond[2]
  stats_df_new$coefs2[[i]] <- fixef(new_mods[[i]])$cond[3]
}

#### TODO - Export table ####
# plot table for now
stats_df_new %>%
  mutate(effects2 = case_when(Taxon == "Anthocoridae" ~ "None",
                              Taxon != "Anthocoridae" ~  effects2, )) %>%
  mutate(coefs2 = case_when(Taxon == "Anthocoridae" ~ 0,
                              Taxon != "Anthocoridae" ~  coefs2, )) %>%
  group_by(Taxon) %>%
  tab_df(title = "Top predator models - vegdata")
# note: these models are better than the corresponding "no veg" models

# Extra ladybug model ####
# try a binomial cocc mod for fall (low coc density in fall)
# source: build models for each distweight and dredge
source("coccinellidae_binomial.R", echo = TRUE)

# # check mod table
# all_bin_mods %>% View

# # review_mod tables and top mods
# all_bin_mods %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>% View


# Aphid models ####
## Spring
# AllAph spring
r2(get.models(tab_nb_allaph_sp_scaled, 1)[[1]]) # good fit 0.8
plot(simulateResiduals(get.models(tab_nb_allaph_sp_scaled, 1)[[1]])) # ok
tab_nb_allaph_sp_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View
# wateringMethod in all top mods, which are close in deltas w/low weights
# -ag_sig1, +alfalfa_sig1, +Cocc

# Acrythosiphon spring
r2(get.models(tab_nb_acy_sp_scaled, 1)[[1]]) # good fit 0.8
plot(simulateResiduals(get.models(tab_nb_acy_sp_scaled, 1)[[1]])) # ok
tab_nb_acy_sp_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View
# top mods generally have different landcover factors, but all at _sig1 scale
# watering method still in all top mods
# model average?

# nonacy spring
r2(get.models(tab_nb_nonacy_sp_scaled, 1)[[1]]) # average fit 0.56
plot(simulateResiduals(get.models(tab_nb_nonacy_sp_scaled, 1)[[1]])) # great
tab_nb_nonacy_sp_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View
# top mods all include -dirt_no

# Spring Summary
## SEM should include coccinellidae, wateringmethod for sure.
## maybe also include ag_sig1, alfalfa_sig1

### Fall
# AllAph fall
r2(get.models(tab_nb_allaph_fa_scaled, 1)[[1]]) # average fit 0.69
plot(simulateResiduals(get.models(tab_nb_allaph_fa_scaled, 1)[[1]])) # ok
tab_nb_allaph_fa_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View

# clear top mod with +ich, -impermeable_sig4, +naturalArid_sig4

# Acyrthosiphon fall
r2(get.models(tab_nb_acy_fa_scaled, 1)[[1]]) # pretty good fit 0.7
plot(simulateResiduals(get.models(tab_nb_acy_fa_scaled, 1)[[1]])) # great
tab_nb_acy_fa_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View
# +ich in all mods. -geo in top 3 mods. +naturalArid across scales!

# Nonacy fall
r2(get.models(tab_nb_nonacy_fa_scaled, 1)[[1]]) # pretty good fit 0.65
plot(simulateResiduals(get.models(tab_nb_nonacy_fa_scaled, 1)[[1]])) # weird?
tab_nb_nonacy_fa_scaled %>%
  tibble %>%
  slice(1:5) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  View
# +ich in all mods. landcover effects varied in top mods,
# BUT top mod is way better than #2. Includes -impermeable_sig2, +natArid_sig2

# Fall Summary
# strong correlations with ichneumonoidea.
# natural arid benefits aphid abundance. _sig4 best scale for AllAph.

# Review aphid models ####
# spring - #1 mod is inappropriate because it combines wateringmethod and
# water_sig1
sp_best <- get.models(tab_nb_allaph_sp_scaled, 2)[[1]]
fa_best <- get.models(tab_nb_allaph_fa_scaled, 1)[[1]]
summary(sp_best)
summary(fa_best)

## make aphid modstats table
aph_mods <- list("allaphFA" = sp_best, "araFA" = fa_best)
# build empty tibble to hold stats
stats_df_aph <- tibble(Taxon = c("AllAph", "AllAph"),
                      Season = c("Spring", "Fall"),
                      MarginalR2 = c(0),
                      ConditionalR2 = c(0),
                      effects1 = c("none"),
                      effects2 = c("none"),
                      effects3 = c("none"),
                      coefs1 = c(0),
                      coefs2 = c(0),
                      coefs3 = c(0))
# fill tibble with stats
for (i in seq_along(aph_mods)){
  stats_df_aph$MarginalR2[[i]] <- r2(aph_mods[[i]])[[2]]
  stats_df_aph$ConditionalR2[[i]] <- r2(aph_mods[[i]])[[1]]
  stats_df_aph$effects1[[i]] <- names(aph_mods[[i]]$frame)[2]
  stats_df_aph$effects2[[i]] <- names(aph_mods[[i]]$frame)[3]
  stats_df_aph$effects3[[i]] <- names(aph_mods[[i]]$frame)[4]
  stats_df_aph$coefs1[[i]] <- fixef(aph_mods[[i]])$cond[2]
  stats_df_aph$coefs2[[i]] <- fixef(aph_mods[[i]])$cond[3]
  stats_df_aph$coefs3[[i]] <- fixef(aph_mods[[i]])$cond[4]
}

#### TODO - Export table ####
# plot table for now
stats_df_aph %>%
  tab_df

# Aphid models with vegdata ####
# recall vegdata from predator modeling
df_sp_vd
df_fa_vd
# add vegdata to top mods
sp_best
sp_best_veg <- glmmTMB(AllAph ~ Treatment + ag_sig1 + wateringMethod +
                         log(Coccinellidae + 1) + shan + rich + totalCover +
                         (1 | Site:Field),
                       data = df_sp_vd,
                       family = "nbinom2",
                       na.action = "na.fail")

veg_dredge_sp <- dredge(sp_best_veg,
                        fixed = "cond(Treatment)",
                        m.lim = c(0, 4),
                        trace = 2)

fa_best
fa_best_veg <- glmmTMB(AllAph ~ Treatment + impermeable_sig4 +
                         naturalArid_sig4 + log(Ichneumonoidea + 1) + shan +
                         rich + totalCover + (1 | Site:Field),
                    data = df_fa_vd,
                    family = "nbinom2",
                    na.action = "na.fail")

veg_dredge_fa <- dredge(fa_best_veg,
                        fixed = "cond(Treatment)",
                        m.lim = c(0, 4),
                        trace = 2)

# review model selection tables and best models
## Spring
veg_dredge_sp %>% View
# new best mod!! +richness!!
sp_best_veg <- get.models(veg_dredge_sp, 1)[[1]]
summary(sp_best_veg)
r2(sp_best_veg)
plot(simulateResiduals(sp_best_veg))
plot(fitted(sp_best_veg), residuals(sp_best_veg, type = "pearson"))
plot(df_sp_vd$AllAph, fitted(sp_best_veg))
abline(0, 1)

plot(allEffects(sp_best_veg, resid = TRUE))

## Fall
veg_dredge_fa %>% View
# new best mod!! -shan!!
fa_best_veg <- get.models(veg_dredge_fa, 1)[[1]]
summary(fa_best_veg)
r2(fa_best_veg)
plot(simulateResiduals(fa_best_veg))
plot(fitted(fa_best_veg), residuals(fa_best_veg, type = "pearson"))
plot(df_fa_vd$AllAph, fitted(fa_best_veg))
abline(0, 1)

plot(allEffects(fa_best_veg, resid = TRUE))



# Sham attraction ####
# source: identify predators that are attracted to sham treatments via
# mixed-effects models of the *difference* between sham and control plots.
source("sham_attraction.R", echo = TRUE)

## Examine models
diff_stats_sp # Ara, Coc *; Anth ***
diff_stats_fa # none significant

## Examine Ichneumonidae - generalized linear model
summary(d_ich_mod)
#try nb mod for ich in fall
nb_d_ich <- glmmTMB(Ichneumonoidea ~ Treatment + (1 | Site:Field),
                     family = "nbinom2",
                     data = subplot_data_raw %>% filter(Season == "Fall",
                                                      Treatment != "Pre-"))

summary(nb_d_ich)
exp(0.8307) # exponentiation of sham effect estimate
# this is better I think. could do pois or nb model. Overdispersion is there,
# but minimal.



### Plot sham effects
#### TODO - come up with a better plot ####
# df of asterisks
sigs1_sp <- tibble(Taxon = c("Anthocoridae"),
               Difference = rep(16, 1),
               Season = c("Spring"))
sigs2_sp <- tibble(Taxon = c("Arachnida", "Coccinellidae"),
                Difference = rep(16, 2),
                Season = c("Spring"))


plot_w_diff %>% pivot_longer(contains("diff"),
                           names_to = "Taxon",
                           values_to = "Difference") %>%
  mutate(Taxon = str_sub(Taxon, 5)) %>%
  filter(Taxon %in% c("Anthocoridae",
                      "Arachnida",
                      "Coccinellidae",
                      "Geocoris")) %>%
  ggplot(aes(Taxon, Difference)) +
  geom_boxplot() +
  facet_grid(. ~ Season) +
  theme_classic() +
  geom_text(data = sigs1_sp, label = "***") +
  geom_text(data = sigs2_sp, label = "*") +
  geom_hline(yintercept = 0, color = "red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


# the above plot is ugly. make a better one.
plot_w_diff %>% pivot_longer(contains("diff"),
                           names_to = "Taxon",
                           values_to = "Difference") %>%
  mutate(Taxon = str_sub(Taxon, 5)) %>%
  filter(Taxon %in% c("Anthocoridae",
                      "Arachnida",
                      "Coccinellidae",
                      "Geocoris")) %>%
  filter(Season == "Fall") %>%
  ggplot(aes(y = Taxon, x = Difference)) +
  geom_density_ridges(scale = 1.5) +
  # facet_wrap(~Season) +
  geom_vline(xintercept = 0, color = "red") +
  theme_classic()
# this isn"t much better!

# Top-down effects ####
# source: this portion moved to a separate script for clarity
source("top_down.R", echo = TRUE)

# SEM-LAVAAN ####
### after talking to beth...####



# 1. prepare data
library(sjmisc) # to_dummy

lavaan_df_prep <- subplot_data_raw %>%
  to_dummy(wateringMethod, suffix = "label") %>%
  bind_cols(subplot_data_raw)
lavaan_df <- lavaan_df_prep %>%
  to_dummy(Treatment, suffix = "label") %>%
  bind_cols(lavaan_df_prep) %>%
  left_join(diffData_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-",
         diffCoccinellidae > 0) %>%
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) # %>%
  # mutate(diffCoccinellidae =
  #          case_when(Treatment_Sham == 1 ~ diffCoccinellidae,
  #                    Treatment_Sham == 0 ~ 0))


# 2. review aphid ladybug relationship
## from prior results....
cocc_eff <- glmmTMB(AllAph ~ Treatment + log(Coccinellidae + 1) +
                      (1 | Site:Field),
                    data = lavaan_df,
                    family = "nbinom2")
summary(cocc_eff) # -TreatmentSham**
plot(allEffects(cocc_eff))
plot(simulateResiduals(cocc_eff))


# 3. build SEM
## try lavaan ####
library(lavaan)
library(lavaanPlot)

mod_spec <- "
  # direct effects
    AllAph ~ lcONE*wateringMethod_Flooding + lcTWO*ag_sig1 + b*diffCoccinellidae
    Coccinellidae ~ lcTHREE*weedy_sig1 + lcFOUR*dirt_sig1
    diffCoccinellidae ~ c*Coccinellidae + d*Treatment_Sham
  # mediator
    # ?
  # indirect effects
    # Coccinellidae > diffCoccinellidae > AllAph
    cb := c*b
    # Treatment > diffCoccinelliade > AllAph
    db := d*b
  # total effects
    # Treatment > AllAph
    # tTreat := (d*b)
    # Coccinellidae > AllAph
    # ?
  # covariance ## ERROR - model not identifiable, maybe need to fix a param
  # AllAph ~~ Coccinellidae
  # constraints
    d == 1
    lcTHREE == 1
"

# test
mod_spec <- "
Coccinellidae ~ b*weedy_sig1 + dirt_sig1
diffCoccinellidae ~ Coccinellidae + a*Treatment_Sham
AllAph ~ diffCoccinellidae + wateringMethod_Flooding + ag_sig1
# constraint
a == 1
b == 1
"
mod_fit <- sem(mod_spec, data = lavaan_df)
summary(mod_fit)
lavaanPlot(model = mod_fit, coefs = TRUE, covs = FALSE, stand = FALSE)

### TODO - understand scale of coeffs. Check constraints ####
