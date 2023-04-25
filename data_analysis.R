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
"#548235" # Acyrthosiphon aphids
"#3B3838" # Non-Acyrthosiphon aphids
"#76db91" # Spring
"#9e3c21" # Fall:
"Spectral" # Predators (alphabe)
"Set1" # Sites

# Landcover: Use palette derived from actual colors?
# trainingPalette asset from earth engine:
"#173118" # alfalfa
"#7f736d" # naturalArid
"#93837c" # weedy
"#373f2a" # riparian
"#696153" # Ag group midpoint
  "#a38878" # rye
  "#2e3a2d" # onion
"#A8A4A5" # Bare group midpoint (first two)
  "#b6b2b0" # bare
  "#999599" # dirtRoad
  "#6f5a4f" # dairy
"#89939F" # Impermeable group midpoint
  "#50514c" # paved
  "#c2d4f1" # structures
NA # Surface water

# Contrast too low, make curstom palette

lc_palette <- c(
  "#156b07", # alfalfa
  "#cca266", # naturalArid
  "#a5cc66", # weedy
  "#192280", # riparian
  "#cc6655", # Ag group
  "#801f19", # Bare group
  "grey40", # Impermeable
  "#9987f1" # log Aphid density
)
# ensure factors are ordered correctly
lc_fct_ord <- c()

# distance
dist_palette <- c(
  "#fafa6e",
  "#86d780",
  "#23aa8f",
  "#007882",
  "#2a4858"
)

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
                                "Geocoris", "Ichneumonoidea"), 2),
                  Season = c(rep("Spring", 5), rep("Fall", 5)),
                  MarginalR2 = c(0),
                  ConditionalR2 = c(0),
                  effects1 = c("none"),
                  effects2 = c("none"),
                  coefs1 = c(0),
                  coefs2 = c(0),
                  coefs1min = c(0),
                  coefs1max = c(0),
                  coefs2min = c(0),
                  coefs2max = c(0),
                  coefsse1 = c(0),
                  coefsse2 = c(0),
                  )

# fill tibble with stats
for (i in seq_along(best_mod_list)){
  stats_df$MarginalR2[[i]] <- r2(best_mod_list[[i]])[[2]]
  stats_df$ConditionalR2[[i]] <- r2(best_mod_list[[i]])[[1]]
  stats_df$effects1[[i]] <- names(best_mod_list[[i]]$frame)[2]
  stats_df$effects2[[i]] <- names(best_mod_list[[i]]$frame)[3]
  stats_df$coefs1[[i]] <- fixef(best_mod_list[[i]])$cond[2]
  stats_df$coefs2[[i]] <- fixef(best_mod_list[[i]])$cond[3]
  stats_df$coefs1min[[i]] <- confint(best_mod_list[[i]])[2,1]
  stats_df$coefs2min[[i]] <- confint(best_mod_list[[i]])[3,1]
  stats_df$coefs1max[[i]] <- confint(best_mod_list[[i]])[2,2]
  stats_df$coefs2max[[i]] <- confint(best_mod_list[[i]])[3,2]
  stats_df$coefsse1[[i]] <- summary(best_mod_list[[i]])$coefficients$cond[2,2]
  stats_df$coefsse2[[i]] <- summary(best_mod_list[[i]])$coefficients$cond[3,2]

}


stats_df %>%
  group_by(Season) %>%
  tab_df(title = "Top predator models (no vegetation data included)")
## BOOKMARK ####
# maybe build the figure in notebook p. 28 based on the table above.
# no spring pred mods need vegdata, so you can rock with this.
# fall mods that need vegdata:
# Anthocoridae, Arachnida

# reshape stats table
figDat <- stats_df %>%
  rename(coefsmin1 = coefs1min,
         coefsmax1 = coefs1max,
         coefsmin2 = coefs2min,
         coefsmax2 = coefs2max) %>%
  pivot_longer(effects1:coefsse2,
               names_to = c(".value", "effectRank"),
               names_pattern = "(.*)(.s*)",
               values_to = c("var1, var2, var3, var4, var5, var6")) %>%
  mutate(effects = case_when(effects == "log_AllAph" ~ "logAllAph_NA",
                             effects == "wateringMethod" ~ "wateringMethod_NA",
                             effects != c("log_AllAph",
                                          "wateringMethod") ~ effects)) %>%
  separate(effects, into = c("Effect", "Distance"), sep = "_") %>%
  mutate(Taxon = fct_rev(as_factor(Taxon))) %>%
  mutate(Distance = factor(Distance, levels = c("NA",
                                                "sig1",
                                                "sig2",
                                                "sig3",
                                                "sig4",
                                                "sig5",
                                                "const",
                                                "no")))


# test ggplot
figDat %>%
  filter(Season == "Spring") %>%
  ggplot(aes(x = Distance, y = Taxon, group = effectRank)) +
  geom_point(aes(color = Effect,
                 size = abs(coefs),
                 shape = as_factor(sign(coefs))),
             position = position_dodge(-0.5),
             stroke = 2) +
  scale_size(range = c(5, 24), name = "Effect size") +
  scale_x_discrete(drop = FALSE) +
  scale_shape_manual(name = "Effect size",
                     values = c(1, 16)) +
  scale_color_brewer(palette = "Dark2")

# CAPTION (v1): Effects of explanatory variables (including land cover features
# and log aphid density) and their distance weighting (x axis) on the density
# of key predators (y axis). Only the explanatory variables from the best model
# of each predator (by lowest AICc) are shown. The identities of the
# explanatory variables are represented by color, and the points are scaled to
# the effect size. Open circles represent negative effect sizes.



# test ggplot v2
figDat %>%
  filter(Season == "Spring") %>%
  mutate(Effect = reorder(Effect, desc(coefs))) %>%
  ggplot(aes(x = Taxon, y = coefs, group = effectRank, color = Distance)) +
  geom_point(stroke = 2.5) +
  geom_linerange(aes(ymin = coefsmin, ymax = coefsmax), size = 1.5) +
  facet_grid(~Effect, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("NA" = "grey50",
                                "sig1" = "#003366",
                                "sig2" = "#174978",
                                "sig3" = "#2f5f8a",
                                "sig4" = "#46769b",
                                "sig5" = "#5e8cad",
                                "const" = "purple",
                                "no" = "violet"
                                ), drop = FALSE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("Effect size")

# CAPTION (v2): Effects of explanatory variables (including land cover features
# and log aphid density) on the density of key predators. Only explanatory
# variables from the best models (lowest AICc) are shown. The distance
# weighting of the effect (where applicable) is denoted by color. Explanatory
# variables are ordered (left to right) by mean effect size across all taxa.
# Error bars represent 95% confidence intervals.

# test ggplot v3
figDat %>%
  filter(Season == "Spring") %>%
  mutate(Effect = reorder(Effect, desc(coefs))) %>%
  ggplot(aes(x = Effect, y = coefs, group = effectRank, color = Distance)) +
  geom_point(stroke = 2.5) +
  geom_linerange(aes(ymin = coefsmin, ymax = coefsmax), size = 1.5) +
  facet_grid(~Taxon, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("NA" = "grey50",
                                "sig1" = "#003366",
                                "sig2" = "#174978",
                                "sig3" = "#2f5f8a",
                                "sig4" = "#46769b",
                                "sig5" = "#5e8cad",
                                "const" = "purple",
                                "no" = "violet"
  ), drop = FALSE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("Effect size")

# CAPTION (v3): Effects of explanatory variables (including land cover features
# and log aphid density) on the density of key predators. Only explanatory
# variables from the best models (lowest AICc) are shown. The distance
# weighting of the effect (where applicable) is denoted by color.
# Error bars represent 95% confidence intervals.

# # test ggplot v4 - no go
# figDat %>%
#   filter(Season == "Spring") %>%
#   ggplot(aes(x = Distance, y = Taxon, group = effectRank)) +
#   geom_point(aes(color = Effect, # estimate
#                  size = (abs(coefs))*sqrt(2),
#                  shape = as_factor(sign(coefs))),
#                  stroke = 1,
#              position = position_dodge(-0.5),
#              alpha = 1) +
#   geom_point(aes(color = Effect, # uncertainty
#                  size = (abs(coefs)-3*estRange)*sqrt(2),
#                  stroke = 6*estRange),
#              shape = 21,
#              position = position_dodge(-0.5),
#              alpha = 0.5) +
#   scale_x_discrete(drop = FALSE) +
#   scale_shape_manual(name = "Effect size",
#                      values = c(13, 1)) +
#   scale_color_brewer(palette = "Dark2")

# test ggplot v5
install.packages("devtools")
library(devtools)
devtools::install_github("wilkelab/ungeviz")
library(ungeviz)

## need to calculate "moe" stat for ungeviz::stat_confidence_density

# both
figDat %>%
  mutate(moe95 = coefsse*1.96,
         Distance = recode(Distance, sig1 = "75", sig5 = "650", "NA" = "0",
                           no = "1000+",
                           sig3 = "350"),
         Effect = recode_factor(Effect,
                                alfalfa = "Alfalfa",
                                natArid = "Desert shrub",
                                weedy = "Weedy cover",
                                riparian = "Riparian",
                                ag = "Non-alfalfa agriculture",
                                dirt = "Bare soil",
                                impermeable = "Impermeable surfaces",
                                logAllAph = "log(Aphid density)",
                                wateringMethod = "Flood irrigation",
                                .ordered = TRUE # ensure factor order to match palette
         ),
         Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
         Taxon = recode_factor(Taxon,
                               Anthocoridae = "Anthocor.",
                               Arachnida = "Arachnida",
                               Coccinellidae = "Coccinell.",
                               Geocoris = "Geocoris",
                               Ichneumonoidea = "Ichneum."),
         roundM = round(MarginalR2, 2),
         quo = "{R^2}[M]~'='",
         cat = paste0(quo, "~'", roundM, "'"),
         expr2 = list(bquote("{R^2}[M]~'='"~.(roundM))),
         expr = paste0(round(MarginalR2,2))) %>%
  # pull(cat)
  ggplot(aes(y = Distance,
             x = coefs,
             group = effectRank,
             fill = Effect,
             # color = Effect,
             moe = moe95)) +
  facet_grid(Taxon~Season) +
  stat_confidence_density(height = 0.9,
                          position = position_dodge(0.8),
                          show.legend = TRUE) +
  # coord_flip() +
  # scale_y_discrete(drop = FALSE) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(xmin = coefs+coefsse,
                    xmax = coefs-coefsse),
                position = position_dodge(0.8),
                width = 0.2
  ) +
  xlim(c(-3.4,3)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red",
             alpha = 0.5) +
  geom_hline(yintercept = (1:4)+0.5,
             color = "grey30") +
  xlab("Standardized effect size") +
  ylab("Spatial scale of effect (m)") +
  scale_fill_manual(name = "Explanatory variable",
                    values = c(
                      "#156b07", # alfalfa
                      # "#cca266", # naturalArid
                      "#a5cc66", # weedy
                      # "#192280", # riparian
                      "#ebad02", # Ag group
                      "#801f19", # Bare group
                      "grey40", # Impermeable
                      "#9987f1", # log Aphid density
                      "#0000ff" # Watering Method
                    )) +
  theme_grey(base_size = 10) +
  theme(#legend.position = c(0.1, 0.1),
        legend.background = element_rect(linetype = 1, color = NA),
        panel.background = element_rect(fill = NA, color = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.1, color = "grey80"),
        panel.grid.minor.x = element_line(size = 0.1, color = "grey90"),
        strip.background.x = element_rect(fill = "NA", color = "NA"),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin =  margin(r = 0.2, l = -40, t = 0)) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_label(aes(x = -3,
                 y = 2.9,
                 label = cat),
             inherit.aes = F,
             parse = T,
             size = 3,
             hjust = 0)

# both - no distance
figDat %>%
  mutate(moe95 = coefsse*1.96,
         Distance = recode(Distance, sig1 = "75", sig5 = "650", "NA" = "0",
                           no = "1000+",
                           sig3 = "350"),
         Effect = recode_factor(Effect,
                                alfalfa = "Alfalfa",
                                natArid = "Desert shrub",
                                weedy = "Weedy cover",
                                riparian = "Riparian",
                                ag = "Non-alfalfa agriculture",
                                dirt = "Bare soil",
                                impermeable = "Impermeable surfaces",
                                logAllAph = "log(Aphid density)",
                                wateringMethod = "Flood irrigation",
                                .ordered = TRUE # ensure factor order to match palette
         ),
         Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
         Taxon = recode_factor(Taxon,
                               Anthocoridae = "Anthocor.",
                               Arachnida = "Arachnida",
                               Coccinellidae = "Coccinell.",
                               Geocoris = "Geocoris",
                               Ichneumonoidea = "Ichneum."),
         roundM = round(MarginalR2, 2),
         quo = "{R^2}[M]~'='",
         cat = paste0(quo, "~'", roundM, "'"),
         expr2 = list(bquote("{R^2}[M]~'='"~.(roundM))),
         expr = paste0(round(MarginalR2,2))) %>%
  # pull(cat)
  ggplot(aes(y = Taxon,
             x = coefs,
             fill = Effect,
             # color = Effect,
             moe = moe95)) +
  facet_grid(~Season) +
  stat_confidence_density(height = 0.9,
                          position = position_dodge(0.8),
                          show.legend = TRUE) +
  coord_flip() +
  # scale_y_discrete(drop = FALSE) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(xmin = coefs+coefsse,
                    xmax = coefs-coefsse),
                position = position_dodge(0.8),
                width = 0.2
  ) +
  xlim(c(-3.4,3)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red",
             alpha = 0.5) +
  geom_hline(yintercept = (1:4)+0.5,
             color = "grey70") +
  xlab("Standardized effect size") +
  ylab("Predator Taxon") +
  scale_fill_manual(name = "Explanatory variable",
                    values = c(
                      "#156b07", # alfalfa
                      # "#cca266", # naturalArid
                      "#a5cc66", # weedy
                      # "#192280", # riparian
                      "#ebad02", # Ag group
                      "#801f19", # Bare group
                      "grey40", # Impermeable
                      "#9987f1", # log Aphid density
                      "#0000ff" # Watering Method
                    )) +
  theme_grey(base_size = 10) +
  theme(#legend.position = c(0.1, 0.1),
    legend.background = element_rect(linetype = 1, color = NA),
    panel.background = element_rect(fill = NA, color = "black"),
    plot.background = element_rect(fill = "white"),
    # panel.grid.major.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background.x = element_rect(fill = "NA", color = "NA"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box.margin =  margin(r = 0.2, l = -40, t = 0)) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_text(aes(x = 3,
                 y = Taxon,
                 label = cat),
             inherit.aes = F,
             parse = T,
             size = 2.5)


ggsave("predmods.png", width = 6.5, height = 4, units = "in")

## additional plot for distance

# define raster grob with horizontal gradient
g <- rasterGrob(matrix(dist_palette, nrow = 1),
                width=unit(1,"npc"), height = unit(1,"npc"),
                interpolate = TRUE)
figDat %>%
  mutate(Distance = recode(Distance, sig1 = "75", sig5 = "650", "NA" = "0",
                           no = "1000",
                           sig3 = "350"),
         Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
         Taxon = recode_factor(Taxon,
                               Anthocoridae = "Anthocoridae",
                               Arachnida = "Arachnida",
                               Coccinellidae = "Coccinellidae",
                               Geocoris = "Geocoris",
                               Ichneumonoidea = "Ichneumonoidea"),
         combo = paste0(Season, Taxon, Distance),
         dnum = as.numeric(as.character(Distance))) %>%
  distinct(combo, .keep_all = TRUE) %>%
  filter(dnum > 0) %>%
  mutate(hjust = c(0,-0.1,-0.2,-0.2,-0.3,
                   0,0,1,1,-0.1),
         y = c(9,7,5,9,3,9,9,9,7,7)) %>%
  select(Season, Taxon, dnum, hjust, y) -> fd2
  ggplot(aes(dnum)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-10, 1040),
                     breaks = c(75, 350, 650, 1000)) +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 6)) +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_vline(aes(xintercept=dnum), color="black", size=1) +
  facet_wrap(~Season, ncol = 1) +
  theme_grey(base_size = 14) +
  geom_label(aes(y = y, label = Taxon, hjust = hjust)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background.x = element_rect(fill = "NA", color = "NA")) +
  xlab("Scale of land cover effects, m")

Distance <- seq(1,1000)
Vaggressive <- c((1.0 / (1.0 + 2.718282 ^ ((3.5 * Distance / 50) - 4.5))) * 1000, rep(0, 100))
Aggressive <- c((1.0 / (1.0 + 2.718282 ^ ((2.0 * Distance / 50) - 6.0))) * 1000, rep(0, 100))
Moderate <- c((1.0 / (1.0 + 2.718282 ^ ((1.25 * Distance / 50) - 8.0))) * 1000, rep(0, 100))
Weak <- c((1.0 / (1.0 + 2.718282 ^ ((1.0 * Distance / 50) - 10.0))) * 1000, rep(0, 100))
Vweak <- c((1.0 / (1.0 + 2.718282 ^ ((0.75 * Distance / 50) - 11.0))) * 1000, rep(0, 100))
Constant <- c((abs(Distance-1000)), rep(0, 100))
None <- c(rep(1000, 1000), rep(0, 100))
# update dist to match length of other cols
Distance <- seq(1, 1100)
df <- data.frame(Distance, Vaggressive, Aggressive, Moderate, Weak, Vweak, Constant, None) %>%
  pivot_longer(-Distance, names_to = "Decay function", values_to = "value") %>%
  mutate(`Decay function` = fct_relevel(`Decay function`, "Vaggressive", "Aggressive", "Moderate", "Weak", "Vweak", "Constant", "None")) %>%
  mutate(`Decay function` = fct_recode(`Decay function`, `Very aggressive` = "Vaggressive", `Very weak` = "Vweak"))

# distance
dist_palette <- c(
  "#fafa6e",
  "#86d780",
  "#23aa8f",
  "#007882",
  "#2a4858",
  "grey50",
  "grey80"
)
figDat %>%
  mutate(Distance = recode(Distance, sig1 = "75", sig5 = "650", "NA" = "0",
                           no = "1000",
                           sig3 = "350"),
         Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
         Taxon = recode_factor(Taxon,
                               Anthocoridae = "Anthocoridae",
                               Arachnida = "Arachnida",
                               Coccinellidae = "Coccinellidae",
                               Geocoris = "Geocoris",
                               Ichneumonoidea = "Ichneumonoidea"),
         combo = paste0(Season, Taxon, Distance),
         dnum = as.numeric(as.character(Distance))) %>%
  distinct(combo, .keep_all = TRUE) %>%
  filter(dnum > 0) %>%
  mutate(hjust = c(0.03,-0.02,-0.07,-0.13,-0.09,
                   0.03,-0.8,1.02,1.02,-0.02),
         y = c(9,7,5,8,3,
               9,4,9,7,7)) %>%
  select(Season, Taxon, dnum, hjust, y) -> fd2
fd2 %>%
  mutate(col = case_when(dnum == 75 ~ "Very aggressive",
                         dnum == 350 ~ "Moderate",
                         dnum == 650 ~ "Very weak",
                         dnum == 1000 ~ "None"),
         col2 = as.factor(col)) %>%
  mutate(col2 = fct_expand(col2, "Aggressive", "Weak", "Constant", "None")) %>%
  mutate(col2 = fct_relevel(col2, "Very aggressive", "Aggressive", "Moderate",
                            "Weak", "Very weak", "Constant", "None"))-> df3
ggplot(df, aes(Distance, value, color = `Decay function`)) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_manual(values = dist_palette) +
  scale_fill_manual(values = dist_palette, drop = F, guide = "none") +
  labs(y = "Weighting value") +
  # geom_vline(data = fd2, aes(xintercept = dnum)) +
  geom_label(data = df3, aes(x = dnum,
                             y = y*100,
                             hjust = hjust,
                             label = Taxon,
                             fill = col2),
             inherit.aes = FALSE,
             alpha = 0.8) +
  facet_wrap(~Season, ncol = 1) +
  theme(strip.background = element_blank())


## next idea
# split info above into one "heatmap" figure" (taxon*predictor, heat = effect size)
# and one distance*taxon figure?




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
# 1. prepare data
## binary watering method var
lavaan_df_prep <- subplot_data_raw %>%
  to_dummy(wateringMethod, suffix = "label") %>%
  bind_cols(subplot_data_raw)
## binary site vars
lavaan_df_prep2 <- lavaan_df_prep %>%
  to_dummy(Site, suffix = "label") %>%
  bind_cols(lavaan_df_prep)
## binary treatment vars, bind other data

## one version only - you get lost if you make multiple dfs!
lavaan_df1 <- lavaan_df_prep2 %>%
  to_dummy(Treatment, suffix = "label") %>%
  bind_cols(lavaan_df_prep2) %>%
  left_join(diffData_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-") %>%
         # diffCoccinellidae > 0) %>%
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) %>%
  select(logAllAph, logCoccinellidae, Coccinellidae, # dump other cols
         diffCoccinellidae,
         wateringMethod_Flooding,
         ag_sig1, weedy_sig1, dirt_sig1,
         AllAph, Treatment,
         Treatment_Sham, Site, Field) %>%
  mutate(diffCoccinellidae =
           case_when(Treatment_Sham == 1 ~ diffCoccinellidae,
                     Treatment_Sham == 0 ~ 0)) %>%
  mutate(across(.cols = -c(Treatment, Site, Field, Treatment_Sham,
                           wateringMethod_Flooding),
                .fns = ~as.vector(scale(.x))))
                  # alternative: ~(. -mean(.)/ sd(.))


## explore data
names(lavaan_df1)
nrow(lavaan_df1)
mean(lavaan_df1 %>% pull(diffCoccinellidae))


# explore dotcharts of all vars
dotchart(lavaan_df1$logAllAph)
dotchart(lavaan_df1$logCoccinellidae)
dotchart(lavaan_df1$diffCoccinellidae)
dotchart(lavaan_df1$logCoccinellidae)
dotchart(lavaan_df1$ag_sig1)
dotchart(lavaan_df1$weedy_sig1)
dotchart(lavaan_df1$dirt_sig1)
dotchart(lavaan_df1$Treatment_Sham)
dotchart(lavaan_df1$wateringMethod_Flooding)
dotchart(subplot_data_raw$weedy_sig1/lala)
dotchart(subplot_data_raw$weedy_sig1)
sd(subplot_data_raw$weedy_sig1)-> lala
# 2. review aphid/ladybug/landcover relationships from glm
## from "top-down" effects model....
cocc_eff <- lmer(logAllAph ~ Treatment + logCoccinellidae +
                      (1 | Site:Field),
                    data = lavaan_df1)
summary(cocc_eff) # !!! Treatment_Sham N.S. unless diffCoccinellidae > 0
plot(allEffects(cocc_eff))
plot(simulateResiduals(cocc_eff))

## from spring allaph model....
summary(get.models(tab_nb_allaph_sp_scaled, 1)[[1]])

## from spring coccinellidae model....
summary(best_mod_list$best.coc.sp)

# 3. build SEM
## try lavaan ####
## 3a. Build model on data subset to estimate sham_trt effect
## subset data (before scaling)
## use data subset to find treatment -> aphid effect
lavaan_df2 <- lavaan_df_prep2 %>%
  to_dummy(Treatment, suffix = "label") %>%
  bind_cols(lavaan_df_prep2) %>%
  left_join(diffData_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-",
         diffCoccinellidae > 0) %>% # filter here!
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) %>%
  select(logAllAph, logCoccinellidae, Coccinellidae, # dump other cols
         diffCoccinellidae,
         wateringMethod_Flooding,
         ag_sig1, weedy_sig1, dirt_sig1,
         AllAph, Treatment,
         Treatment_Sham, Site, Field) %>%
  mutate(diffCoccinellidae =
           case_when(Treatment_Sham == 1 ~ diffCoccinellidae,
                     Treatment_Sham == 0 ~ 0)) %>%
  mutate(across(.cols = -c(Treatment, Site, Field,
                           wateringMethod_Flooding,
                           Treatment_Sham, diffCoccinellidae),
                .fns = ~(./sd(.)))) # scale, don't center
## specify model (with no constraints)
mod_specA <- "
  # direct effects
    logCoccinellidae ~ w*weedy_sig1 + d*dirt_sig1

    logAllAph ~ f*wateringMethod_Flooding + a*ag_sig1 + t*Treatment_Sham

  # mediator
    # ?

  # indirect effects
    # ?
    # multiply coefs here - use :=

  # total effects
    wd := w+(-d) # landcover -> ladybugs
    fa := f+(-a) # landcover -> aphids
    # add coefs here - use :=

  # covariance
   logAllAph ~~ l*logCoccinellidae

  # constraints
   # none yet!
"

## fit model
mod_fitA <- sem(mod_specA, data = lavaan_df2)

summary(mod_fitA)
# Trt_Sham effect is:
# Trtmnt_Shm (t) Est: -0.758    Std.err: 0.133   z: -5.708    P: 0.000

lavaanPlot(model = mod_fitA, coefs = TRUE, covs = F, stand = FALSE)


## 3b. refit to full dataset with Trt_Sham coef "fixed"
mod_specB <- "
  # direct effects
    logCoccinellidae ~ w*weedy_sig1 + d*dirt_sig1

    logAllAph ~ f*wateringMethod_Flooding + a*ag_sig1 + t*Treatment_Sham

  # mediator
    # ?

  # indirect effects
    # ?
    # multiply coefs here - use :=

  # total effects
    wd := w+(-d) # landcover -> ladybugs
    fa := f+(-a) # landcover -> aphids
    # add coefs here - use :=

  # covariance
   logAllAph ~~ l*logCoccinellidae

  # constraints
   t == -0.758
"

# ML is default estimator
mod_fit1B <- sem(mod_specB, data = lavaan_df1)
summary(mod_fit1B)

lavaanPlot(model = mod_fit1B, coefs = TRUE, covs = TRUE, stand = FALSE)


lavResiduals(mod_fit1B)



# exploratory regression
## how does filtering change the sham > diff relationship?
mod1 <- lm(diffCoccinellidae ~ Treatment_Sham + logCoccinellidae,
           data = lavaan_df1)
mod2 <- lm(diffCoccinellidae ~ Treatment_Sham + logCoccinellidae,
           data = lavaan_df2)

summary(mod1)
summary(mod2) # better fit, stronger sham effect.

plot(allEffects(mod1))
plot(allEffects(mod2))

# what is the deal with the diffCoccinellidae



# for negbin SEM - maybe do blavaan, export jags model syntax, edit link
# function, fit w/rjags. nightmarish.

## back to sem?
library(piecewiseSEM)

model <- psem(lm(logCoccinellidae ~ weedy_sig1 + dirt_sig1, lavaan_df2),
              lm(logAllAph ~ wateringMethod_Flooding + ag_sig1 +
                   Treatment_Sham*logCoccinellidae, lavaan_df2))
summary(model)
plot(model)

library(MASS)
# glm version
glm_model <- psem(glm.nb(Coccinellidae ~ weedy_sig1 +
                             dirt_sig1, lavaan_df2),
              glm.nb(AllAph ~ wateringMethod_Flooding + ag_sig1 +
                   Treatment_Sham*Coccinellidae, lavaan_df2))
summary(glm_model)
plot(glm_model)


