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
library(devtools) # for installing packages (ungeviz) from github
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
library(ungeviz) # For "confidence density" ggplots
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

# import/wrangle data ####
# subplot-level data, all vars
# pre- data not /3, can use area offset to correct.
# No transformations have been applied
subplot_data_raw <- read_csv("tidy_data/subplot_data_raw.csv",
                           col_types = "ffffff")

# "Pre-" density has already been /3
# No transformations have been applied
subplot_data <- read_csv("tidy_data/subplot_data.csv",
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
aphids_palette <- c(
  "#548235", # Acyrthosiphon aphids
  "#3B3838" # Non-Acyrthosiphon aphids
)
seasons_palette <- c(
  "#76db91", # Spring
  "#9e3c21" # Fall:
)
"Spectral" # Predators (alphabetical sort)
"Set1" # Sites

# landcover: Use palette derived from actual colors?
lc_palette_experimental <- c(
  "#105405", # alfalfa
  "#aecf7a", # weedy cover
  "#5b3769", # non-alfalfa  ag
  "#801f19", # Bare group
  "grey40", # Impermeable
  "#ff6d00", # log(aphid dens.)
  "#393ebf", # flood irrigation
  "#ff2d55",  # total cover (veg survey)
  "#AAAA20",  # Shannon diversity (land cover)
  "#5856d6" # riparian
  # "#0000FF", # naturalAridPerimeter
  # "#FF0000" # impermeablePerimeter
  # "#cca266" # naturalArid
)

# ensure factors are ordered correctly!
lc_fct_ord <- c()

# distance
dist_palette <- c(
  "#fafa6e", # sig 1
  "#86d780", # sig 2
  "#23aa8f", # sig 3
  "#007882", # sig 4
  "#2a4858", # sig 5
  "grey50",  # const
  "grey80"   # no
)

# Data exploration ####
### SOURCE ####
## Arthropod data summary
source("data_analysis_arthropod_exploration.R", echo = TRUE)

## Landcover data summary
source("data_analysis_classification_summary.R", echo = TRUE)



# Wrangle data for model selection ####
# regular (no rank transform)
# spring

df_sp <- subplot_data_raw %>%
  # spring only
  filter(Season == "Spring") %>%
  # "regular" landcover only
  select(!contains("fix")) %>% #"fix" landcover has manually fixed alfalfa areas
  # log-transform AllAph col
  mutate(log_AllAph = log(AllAph + 1),
         log_Anthocoridae = log(Anthocoridae + 1),
         log_Arachnida = log(Arachnida + 1),
         log_Geocoris = log(Geocoris + 1),
         log_Ichneumonoidea = log(Ichneumonoidea + 1),
         log_Coccinellidae = log(Coccinellidae + 1)) %>%

  # rowwise() %>%
  # # change areascores to proportions
  # mutate(
  #        across(ends_with("_sig1"), ~ .x/rowSums(across(ends_with("_sig1")))),
  #        across(ends_with("_sig2"), ~ .x/rowSums(across(ends_with("_sig2")))),
  #        across(ends_with("_sig3"), ~ .x/rowSums(across(ends_with("_sig3")))),
  #        across(ends_with("_sig4"), ~ .x/rowSums(across(ends_with("_sig4")))),
  #        across(ends_with("_sig5"), ~ .x/rowSums(across(ends_with("_sig5")))),
  #        across(ends_with("_const"), ~ .x/rowSums(across(ends_with("_const")))),
  #        across(ends_with("_no"), ~ .x/rowSums(across(ends_with("_no")))),
  #        ) %>%
  # make diversity cols
  mutate(div_sig1 = diversity(across(ends_with("sig1")), index = "simpson"),
         div_sig2 = diversity(across(ends_with("sig2")), index = "simpson"),
         div_sig3 = diversity(across(ends_with("sig3")), index = "simpson"),
         div_sig4 = diversity(across(ends_with("sig4")), index = "simpson"),
         div_sig5 = diversity(across(ends_with("sig5")), index = "simpson"),
         div_const = diversity(across(ends_with("const")), index = "simpson"),
         div_no = diversity(across(ends_with("no")), index = "simpson"),
         divShan_sig1 = diversity(across(ends_with("sig1")), index = "shannon"),
         divShan_sig2 = diversity(across(ends_with("sig2")), index = "shannon"),
         divShan_sig3 = diversity(across(ends_with("sig3")), index = "shannon"),
         divShan_sig4 = diversity(across(ends_with("sig4")), index = "shannon"),
         divShan_sig5 = diversity(across(ends_with("sig5")), index = "shannon"),
         divShan_const = diversity(across(ends_with("const")),
                                   index = "shannon"),
         divShan_no = diversity(across(ends_with("no")), index = "shannon")
         )  %>%
  # ungroup() %>%
  # center and scale (not needed with rank transform) all landcover + log_AllAph
  # + log_predators
  # hold off on this until later
  mutate(across(.cols = contains(c("_", "shan", "rich", "totalCover", "Perim")),
                .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

## special - fixed only
df_sp_fixed <- subplot_data_raw %>%
  # spring only
  filter(Season == "Spring") %>%
  # "regular" landcover only
  select(!(contains(c("sig1",
                      "sig2",
                      "sig3",
                      "sig4",
                      "sig5",
                      "const",
                      "no")
  ) & !contains("fix")), Ichneumonoidea) %>% #"fix" landcover has manually fixed alfalfa areas
  rename_with(str_replace, pattern = "_fix", replacement = "", .cols=contains("_fix")) %>%
  # log-transform AllAph col
  mutate(log_AllAph = log(AllAph + 1),
         log_Anthocoridae = log(Anthocoridae + 1),
         log_Arachnida = log(Arachnida + 1),
         log_Geocoris = log(Geocoris + 1),
         log_Ichneumonoidea = log(Ichneumonoidea + 1),
         log_Coccinellidae = log(Coccinellidae + 1)) %>%
  # make diversity cols
  mutate(div_sig1 = diversity(across(ends_with("sig1")), index = "simpson"),
         div_sig2 = diversity(across(ends_with("sig2")), index = "simpson"),
         div_sig3 = diversity(across(ends_with("sig3")), index = "simpson"),
         div_sig4 = diversity(across(ends_with("sig4")), index = "simpson"),
         div_sig5 = diversity(across(ends_with("sig5")), index = "simpson"),
         div_const = diversity(across(ends_with("const")), index = "simpson"),
         div_no = diversity(across(ends_with("no")), index = "simpson"),
         divShan_sig1 = diversity(across(ends_with("sig1")), index = "shannon"),
         divShan_sig2 = diversity(across(ends_with("sig2")), index = "shannon"),
         divShan_sig3 = diversity(across(ends_with("sig3")), index = "shannon"),
         divShan_sig4 = diversity(across(ends_with("sig4")), index = "shannon"),
         divShan_sig5 = diversity(across(ends_with("sig5")), index = "shannon"),
         divShan_const = diversity(across(ends_with("const")),
                                   index = "shannon"),
         divShan_no = diversity(across(ends_with("no")), index = "shannon")
  )  %>%
  # ungroup() %>%
  # center and scale (not needed with rank transform) all landcover + log_AllAph
  # + log_predators
  # hold off on this until later
  mutate(across(.cols = contains(c("_", "shan", "rich", "totalCover", "Perim")),
                .fns = ~as.vector(scale(.)))) %>%
  # make area offset
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

df_sp$alfalfa_const == df_sp_fixed$alfalfa_const

dotchart(sort(df_sp$impermeable_sig4), main = "proportion")
# fall
df_fa <- subplot_data_raw %>%
  filter(Season == "Fall") %>%
  select(!contains("fix")) %>%
  mutate(log_AllAph = log(AllAph + 1),
         log_Anthocoridae = log(Anthocoridae + 1),
         log_Arachnida = log(Arachnida + 1),
         log_Geocoris = log(Geocoris + 1),
         log_Ichneumonoidea = log(Ichneumonoidea + 1),
         log_Coccinellidae = log(Coccinellidae + 1)) %>%
  # rowwise() %>%
  # # change areascores to proportions
  # mutate(
  #   across(ends_with("_sig1"), ~ .x/rowSums(across(ends_with("_sig1")))),
  #   across(ends_with("_sig2"), ~ .x/rowSums(across(ends_with("_sig2")))),
  #   across(ends_with("_sig3"), ~ .x/rowSums(across(ends_with("_sig3")))),
  #   across(ends_with("_sig4"), ~ .x/rowSums(across(ends_with("_sig4")))),
  #   across(ends_with("_sig5"), ~ .x/rowSums(across(ends_with("_sig5")))),
  #   across(ends_with("_const"), ~ .x/rowSums(across(ends_with("_const")))),
  #   across(ends_with("_no"), ~ .x/rowSums(across(ends_with("_no")))),
  # ) %>%
  mutate(div_sig1 = diversity(across(ends_with("sig1")), index = "simpson"),
         div_sig2 = diversity(across(ends_with("sig2")), index = "simpson"),
         div_sig3 = diversity(across(ends_with("sig3")), index = "simpson"),
         div_sig4 = diversity(across(ends_with("sig4")), index = "simpson"),
         div_sig5 = diversity(across(ends_with("sig5")), index = "simpson"),
         div_const = diversity(across(ends_with("const")), index = "simpson"),
         div_no = diversity(across(ends_with("no")), index = "simpson"),
         divShan_sig1 = diversity(across(ends_with("sig1")), index = "shannon"),
         divShan_sig2 = diversity(across(ends_with("sig2")), index = "shannon"),
         divShan_sig3 = diversity(across(ends_with("sig3")), index = "shannon"),
         divShan_sig4 = diversity(across(ends_with("sig4")), index = "shannon"),
         divShan_sig5 = diversity(across(ends_with("sig5")), index = "shannon"),
         divShan_const = diversity(across(ends_with("const")),
                                   index = "shannon"),
         divShan_no = diversity(across(ends_with("no")), index = "shannon")
  )  %>%
  # ungroup() %>%
  mutate(across(.cols = contains(c("_", "shan", "rich", "totalCover")),
                .fns = ~as.vector(scale(.)))) %>%
    mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))



# rank transform on landcover vars
# spring
df_sp_rnk <- subplot_data_raw %>%
  filter(Season == "Spring") %>%
  select(!contains("fix")) %>%
  rowwise() %>%
  mutate(div_sig1 = diversity(across(ends_with("sig1")), index = "simpson"),
         div_sig2 = diversity(across(ends_with("sig2")), index = "simpson"),
         div_sig3 = diversity(across(ends_with("sig3")), index = "simpson"),
         div_sig4 = diversity(across(ends_with("sig4")), index = "simpson"),
         div_sig5 = diversity(across(ends_with("sig5")), index = "simpson"),
         div_const = diversity(across(ends_with("const")), index = "simpson"),
         div_no = diversity(across(ends_with("no")), index = "simpson"),
         divShan_sig1 = diversity(across(ends_with("sig1")), index = "shannon"),
         divShan_sig2 = diversity(across(ends_with("sig2")), index = "shannon"),
         divShan_sig3 = diversity(across(ends_with("sig3")), index = "shannon"),
         divShan_sig4 = diversity(across(ends_with("sig4")), index = "shannon"),
         divShan_sig5 = diversity(across(ends_with("sig5")), index = "shannon"),
         divShan_const = diversity(across(ends_with("const")),
                                   index = "shannon"),
         divShan_no = diversity(across(ends_with("no")), index = "shannon")
  )  %>%
  ungroup() %>%
  # rank transform landcover to uniformly distribute values
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~dense_rank(.x))) %>%
  mutate(log_AllAph = log(AllAph + 1)) %>%
  # # center and scale (not needed with rank transform)
  # mutate(across(.cols = contains("_"), # all landcover + log_AllAph
  #               .fns = ~as.vector(scale(.)))) %>%
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

# fall
df_fa_rnk <- subplot_data_raw %>%
  filter(Season == "Fall") %>%
  select(!contains("fix")) %>%
  rowwise() %>%
  mutate(div_sig1 = diversity(across(ends_with("sig1")), index = "simpson"),
         div_sig2 = diversity(across(ends_with("sig2")), index = "simpson"),
         div_sig3 = diversity(across(ends_with("sig3")), index = "simpson"),
         div_sig4 = diversity(across(ends_with("sig4")), index = "simpson"),
         div_sig5 = diversity(across(ends_with("sig5")), index = "simpson"),
         div_const = diversity(across(ends_with("const")), index = "simpson"),
         div_no = diversity(across(ends_with("no")), index = "simpson"),
         divShan_sig1 = diversity(across(ends_with("sig1")), index = "shannon"),
         divShan_sig2 = diversity(across(ends_with("sig2")), index = "shannon"),
         divShan_sig3 = diversity(across(ends_with("sig3")), index = "shannon"),
         divShan_sig4 = diversity(across(ends_with("sig4")), index = "shannon"),
         divShan_sig5 = diversity(across(ends_with("sig5")), index = "shannon"),
         divShan_const = diversity(across(ends_with("const")),
                                   index = "shannon"),
         divShan_no = diversity(across(ends_with("no")), index = "shannon")
  )  %>%
  ungroup() %>%
  mutate(across(.cols = contains("_"), # all landcover + log_AllAph
                .fns = ~dense_rank(.x))) %>%
  mutate(log_AllAph = log(AllAph + 1)) %>%
  mutate(Area = case_when(Treatment == "Pre-" ~ 3,
                          Treatment != "Pre-" ~ 1))

# # example dotcharts - shows distribution of explanatory variables
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp_rnk$AllAph))
dotchart(sort(df_sp$log_AllAph))
dotchart(sort(df_fa$log_AllAph))
dotchart(sort(df_sp$alfalfaPerim))
dotchart(sort(df_sp$naturalAridPerim))
dotchart(sort(df_sp$naturalArid_no))
dotchart(sort(df_sp$impermeable_no))
dotchart(sort(df_sp$impermeablePerim))
dotchart(sort(df_sp$AllAph))
dotchart(sort(df_sp$log_AllAph))
dotchart(sort(df_fa$log_AllAph))


varlist <- c(
  "alfalfa",
  "weedy",
  "weedyWet",
  "ag",
  "dirt",
  "impermeable",
  "water",
  "wet",
  "div",
  "divShan"
)
distlist <- c(
  "_sig1",
  "_sig2",
  "_sig3",
  "_sig4",
  "_sig5",
  "_const",
  "_no"
)
for(i in 1:10){
  for(j in 1:7){
    var <- varlist[[i]]
    dist <- distlist[[j]]
    varDist <- paste0(var, dist)
    dotchart(sort(df_sp[[varDist]]), main = varDist)
  }
}


# Fit models ####

rebuild <- askYesNo("Would you like to rebuild model selection tables?")

if (rebuild == TRUE) {

  ## source external scripts to build model selection tables
  ## set number of cores to be used in parallel processing
  n_cores <- detectCores() - 2
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

# Optional: compare top predator models by viewing model selection tables, ex.
# anth_fams_sp for the top Anthocoridae spring models from each family.


# Make table of top (no veg) predator models ####
# make list of best models
best_mod_list <- list(
  "best.ant.sp" = get.models(nb_scaled$tab_nb_anth_sp_scaled, 1)[[1]],
  "best.ara.sp" = get.models(nb_scaled$tab_nb_ara_sp_scaled, 1)[[1]],
  "best.coc.sp" = get.models(nb_scaled$tab_nb_cocc_sp_scaled, 1)[[1]],
  # "best.coc.sp.fix" = get.models(nb_scaled$tab_nb_cocc_sp_fix_scaled, 1)[[1]],
  "best.geo.sp" = get.models(nb_scaled$tab_nb_geo_sp_scaled, 1)[[1]],
  "best.ich.sp" = get.models(nb_scaled$tab_nb_ich_sp_scaled, 1)[[1]],
  "best.ant.fa" = get.models(nb_scaled$tab_nb_anth_fa_scaled, 1)[[1]],
  "best.ara.fa" = get.models(nb_scaled$tab_nb_ara_fa_scaled, 1)[[1]],
  "best.coc.fa" = get.models(nb_scaled$tab_nb_cocc_fa_scaled, 1)[[1]],
  "best.geo.fa" = get.models(nb_scaled$tab_nb_geo_fa_scaled, 1)[[1]],
  "best.ich.fa" = get.models(nb_scaled$tab_nb_ich_fa_scaled, 1)[[1]]
)

# build empty tibble to hold stats
stats_df <- tibble(Taxon = rep(c("Anthocoridae",
                                 "Arachnida",
                                 "Coccinellidae",
                                 "Geocoris",
                                 "Ichneumonoidea"),
                               2),
                   Season = c(rep("Spring", 5), rep("Fall", 5)),
                   MarginalR2 = c(0),
                   ConditionalR2 = c(0),
                   effects1 = c("none"),
                   effects2 = c("none"),
                   coefs1 = c(0),
                   coefs2 = c(0),
                   coefsmin1 = c(0),
                   coefsmax1 = c(0),
                   coefsmin2 = c(0),
                   coefsmax2 = c(0),
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
  stats_df$coefsmin1[[i]] <- confint(best_mod_list[[i]])[2,1]
  stats_df$coefsmin2[[i]] <- confint(best_mod_list[[i]])[3,1]
  stats_df$coefsmax1[[i]] <- confint(best_mod_list[[i]])[2,2]
  stats_df$coefsmax2[[i]] <- confint(best_mod_list[[i]])[3,2]
  stats_df$coefsse1[[i]] <- summary(best_mod_list[[i]])$coefficients$cond[2,2]
  stats_df$coefsse2[[i]] <- summary(best_mod_list[[i]])$coefficients$cond[3,2]

}

# print table
stats_df %>%
  group_by(Season) %>%
  tab_df(title = "Top predator models (no vegetation data included)")

# Review predator models ####

# best to review these by hand. change input models manually.

# # # optional: review a single mod table
foo <- nb_scaled$tab_nb_cocc_sp_scaled %>%
  tibble %>%
  slice(1:6) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  mutate(data = "RandomForest", .before = everything())
## surprised no perim variables are in here

bar <- tab_nb_cocc_sp_fix_scaled %>%
  tibble %>%
  slice(1:6) %>% # can change how inclusive this is
  select(where(~!all(is.na(.x)))) %>%
  mutate(data = "Manual", .before = everything())

tab_nb_cocc_sp_fix_scaled == nb_scaled$tab_nb_cocc_sp_scaled
importance_tab_fix <- sw(tab_nb_cocc_sp_fix_scaled) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  mutate(names = case_when(names == "cond(Treatment)" ~ "cond(Treatment)_NA",
                           names == "cond(wateringMethod)" ~ "cond(wateringMethod)_NA",
                           names == "cond(log_AllAph)" ~ "cond(logAllAph)_NA",
                           TRUE ~ as.character(names))) %>%
  separate(names, c('class', 'distWeight'), sep = "_") %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `sig1` = 'Very aggressive',
                                       `sig2` = 'Aggressive',
                                       `sig3` = 'Moderately aggressive',
                                       `sig4` = 'Moderate',
                                       `sig5` = 'Slight',
                                       `sig6` = 'Minimal',
                                       `const` = 'Constant',
                                       `no` = 'None',
                                       `NA` = "N/A"))) %>%
  # mutate(distWeight = fct_relevel(distWeight, 'Constant', 'None', after = Inf),
  mutate(class = recode(class,
                        ag = 'Agricultural',
                        alfalfa = 'Alfalfa',
                        dirt = 'Bare soil +\n dirt road',
                        impermeable = 'Impermeable',
                        naturalArid = 'Natural',
                        water = 'Surface\nwater',
                        weedy = 'Weedy',
                        wet = 'Riparian',
                        shan = "Shannon diversity",
                        div = "Simpson diversity"),
         class = fct_relevel(class,
                             "cond(logAllAph)",
                             "cond(Treatment)",
                             "cond(wateringMethod)",
                             .after = "cond(weedy"
                             )) %>%
  rename("sw_fix" = sw, "class_fix" = class, "distWeight_fix" = distWeight)

ggplot(data = importance_tab_fix, aes(x = class_fix, y = distWeight_fix, fill = sw_fix)) +
  geom_tile() +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(angle = 45))+
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       title = paste0('log(Coccinellidae)',
                      ' Variable importance, ',
                      "Spring fixed"))


importance_tab <- sw(nb_scaled$tab_nb_cocc_sp_scaled) %>% #tibble(names = names(.))
  tibble(names = names(.), .name_repair = function(x) gsub('\\.', 'sw', x)) %>%
  arrange(names) %>%
  mutate(names = case_when(names == "cond(Treatment)" ~ "cond(Treatment)_NA",
                           names == "cond(wateringMethod)" ~ "cond(wateringMethod)_NA",
                           names == "cond(log_AllAph)" ~ "cond(logAllAph)_NA",
                           TRUE ~ as.character(names))) %>%
  separate(names, c('class', 'distWeight'), sep = "_") %>%
  mutate(distWeight = as_factor(recode(distWeight,
                                       `sig1` = 'Very aggressive',
                                       `sig2` = 'Aggressive',
                                       `sig3` = 'Moderately aggressive',
                                       `sig4` = 'Moderate',
                                       `sig5` = 'Slight',
                                       `sig6` = 'Minimal',
                                       `const` = 'Constant',
                                       `no` = 'None',
                                       `NA` = "N/A"))) %>%
  # mutate(distWeight = fct_relevel(distWeight, 'Constant', 'None', after = Inf),
  mutate(class = recode(class,
                        ag = 'Agricultural',
                        alfalfa = 'Alfalfa',
                        dirt = 'Bare soil +\n dirt road',
                        impermeable = 'Impermeable',
                        naturalArid = 'Natural',
                        water = 'Surface\nwater',
                        weedy = 'Weedy',
                        wet = 'Riparian',
                        shan = "Shannon diversity",
                        div = "Simpson diversity"),
         class = fct_relevel(class,
                             "cond(logAllAph)",
                             "cond(Treatment)",
                             "cond(wateringMethod)",
                             .after = "cond(weedy"))


ggplot(data = importance_tab, aes(x = class, y = distWeight, fill = sw)) +
  geom_tile() +
  theme(
    axis.text.x=element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(angle = 45))+
  scale_fill_gradient(low="blue", high="red") +
  labs(x = 'Landcover class',
       y = 'Distance weighting algorithm',
       title = paste0('log(Coccinellidae)',
                      ' Variable importance, ',
                      "Spring randomForest"))

# for figure

all_tabs <- cbind(importance_tab_fix, importance_tab) %>%
  mutate(difference = sw_fix - sw,
         distWeight = case_when(
           distWeight == "const)" ~ "Constant",
           distWeight == "no)" ~ "None",
           distWeight == "sig1)" ~ "Proximate",
           distWeight == "sig2)" ~ "Near",
           distWeight == "sig3)" ~ "Intermediate",
           distWeight == "sig4)" ~ "Distant",
           distWeight == "sig5)" ~ "Most Distant",
           distWeight == "N/A" ~ "N/A",
         ),
         class = case_when(
           class == "cond(weedy" ~ "Weedy cover",
           class == "cond(ag" ~ "Non-alfalfa agriculture",
           class == "cond(alfalfa" ~ "Alfalfa",
           class == "cond(dirt" ~ "Bare soil",
           class == "cond(div" ~ "Landcover diversity (Simpson)",
           class == "cond(divShan" ~ "Landcover diversity (Shannon)",
           class == "cond(impermeable" ~ "Impermeable surfaces",
           class == "cond(naturalArid" ~ "Desert shrub",
           class == "cond(water" ~ "Surface water",
           class == "cond(logAllAph)" ~ "log(Aphid density)",
           class == "cond(Treatment)" ~ "Canopy",
           class == "cond(wateringMethod)" ~ "Flood\nirrigation",
         )) %>%
  mutate(distWeight = fct_relevel(distWeight,
    "Proximate",
    "Near",
    "Intermediate",
    "N/A",
    "Distant",
    "Most Distant",
    "Constant",
    "None"
  ),
  class = fct_relevel(class,
                      "Alfalfa",
                      "Bare soil",
                      "Desert shrub",
                      "Flood\nirrigation",
                      "Impermeable surfaces",
                      "Landcover diversity (Simpson)",
                      "Landcover diversity (Shannon)",
                      "log(Aphid density)",
                      "Non-alfalfa agriculture",
                      "Surface water",
                      "Canopy",
                      "Weedy cover"
                      ))

# to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#left side - distance-weighted factors
for_legend <- ggplot(data = all_tabs,
       aes(x = class, y = distWeight, fill = difference)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                      limits = c(-0.3, 0.3)) +
  theme(legend.key.size = unit(25, "pt")) +
  guides(fill = guide_colourbar(title = "Difference in\nvariable importance,\nManual - Random Forest\nalfalfa classification",
                                title.position = "top",
                                direction = "horizontal"))

legend <- g_legend(for_legend)

left <- ggplot(data = all_tabs %>%
                 filter(!distWeight %in% c("Constant", "None", "N/A")),
               aes(x = class, y = distWeight, fill = difference)) +
  geom_tile() +
  theme()+
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       limits = c(-0.3, 0.3)) +
  labs(x = NULL,
       y = 'Distance weighting') +
  theme_grey(base_size = 10) +
  theme(legend.position = "none",
        ### general theme
        panel.background = element_rect(fill = NA, color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 40, hjust = 1),
        strip.background.x = element_rect(fill = "NA", color = "NA"),
        axis.text.y = element_text(angle = 40))


# right side
right <- ggplot(data = all_tabs %>% filter(distWeight == "N/A"),
       aes(x = class, y = distWeight, fill = difference)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       limits = c(-0.3, 0.3)) +
  labs(x = 'Additional predictors\n(Distance weighting N/A)') +
  theme_grey(base_size = 10) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.x=element_text(angle = 45, hjust = 0),
        axis.text.x = element_text(angle = 40, hjust = 1),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background.x = element_rect(fill = "NA", color = "NA")
        )


grid.arrange(left, right, legend,
             layout_matrix = rbind(c(1, 1, 2),
                                   c(1, 1, 2),
                                   c(1, 1, 2),
                                   c(1, 1, 3),
                                   c(1, 1, 3),
                                   c(1, 1, 3),
                                   c(1, 1, 3)))

g <- arrangeGrob(left, right, legend,
             layout_matrix = rbind(c(1, 1, 1, 2),
                                   c(1, 1, 1, 2),
                                   c(1, 1, 1, 2),
                                   c(1, 1, 1, 3),
                                   c(1, 1, 1, 3),
                                   c(1, 1, 1, 3),
                                   c(1, 1, 1, 3)))
ggsave("spring_randomforest_varimportance.pdf",
       g,
       height = 10,
       width = 18,
       units = "cm")




plyr::rbind.fill(foo, bar) %>% tab_df(file = "models.html")
webshot("models.html", "models.pdf")

# choose model to review
review_mod <- get.models(nb_scaled$tab_nb_cocc_sp_scaled, 1)[[1]]

# show summary
summary(review_mod) # no random effect variance. essentially equivalent to
                         # a fixed effects mod. I checked.
# basic effects plots
plot(allEffects(review_mod, residuals = TRUE))

# plot(Effect(c("impermeablePerim"), review_mod, resid = TRUE))
# try fall mod
reviewFigDf <- ggpredict(review_mod, terms = c("weedy_sig1"))
reviewFigDf %>% plot(log.y = T, add.data = T)
# review_mod$frame %>% View

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



fixed_mod <- update(review_mod, data = df_sp_fixed)

summary(review_mod)
summary(fixed_mod)

fixedFigDf <- ggpredict(fixed_mod, terms = c("weedy_sig1"))
reviewFigDf %>% plot(log.y = T, add.data = T)
fixedFigDfAttr <- fixedFigDf %>% attributes()
reviewFigDfAttr <- reviewFigDf %>% attributes()
fixedRaw <- fixedFigDfAttr$rawdata
reviewRaw <- reviewFigDfAttr$rawdata

r2(review_mod)
r2(fixed_mod)

ggplot(fixedRaw, aes(x, response)) +
  geom_point(color = "red") +
  geom_point(data = reviewRaw,
             color = "blue") +
  geom_smooth(aes(x, log(predicted)), data = fixedFigDf, color = "red") +
  geom_smooth(aes(x, log(conf.low)), data = fixedFigDf, color = "red", alpha = 0.5) +
  geom_smooth(aes(x, log(conf.high)), data = fixedFigDf, color = "red", alpha = 0.5) +
  geom_smooth(aes(x, log(predicted)), data = reviewFigDf, color = "blue") +
  geom_smooth(aes(x, log(conf.high)), data = reviewFigDf, color = "blue", alpha = 0.5) +
  geom_smooth(aes(x, log(conf.low)), data = reviewFigDf, color = "blue", alpha = 0.5) +
  geom_text(aes(-1,10, label = "R2m = 0.447\nR2c = NA\nAICc=428.2"), color = "blue") +
  geom_text(aes(-1,9, label = "R2m = 0.380\nR2c = 0.472\nAICc=434.9"), color = "red") +
  labs(title = "marginal effects of weedy_sig1 on ladybug density",
       subtitle = "red = manually classified alfalfa fields\nblue = classified by random forest")

model.sel(review_mod, fixed_mod)

fixedFigDf %>% plot()
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
r2(best.ant.sp.vd) ## new mod has shan instead of impermeable
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
r2(best.ant.fa.vd)
## FORMERLY: total cover looking like the best predictor here
## verdict - OVERTURN. new model is better!
## NOW (w/lc diversity factors): Scale changed from NA to sig3.
## r2 is better in the new model.
## verdict - OVERTURN.

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
## NOW: scale has changed from no to const

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

# plot table
stats_df_new %>%
  group_by(Taxon) %>%
  tab_df(title = "Top predator models - vegdata")
# note: these models are better than the corresponding "no veg" models

# combine info from old no-veg table with this new one
# best_mod_list %>% names

best_mod_list$best.ant.fa <- best.ant.fa.vd
best_mod_list$best.ara.fa <- best.ara.fa.vd

# build empty tibble to hold stats
stats_df_veg <- tibble(Taxon = rep(c("Anthocoridae",
                                     "Arachnida",
                                     "Coccinellidae",
                                     "Geocoris",
                                     "Ichneumonoidea"), 2),
                       Season = c(rep("Spring", 5), rep("Fall", 5)),
                       MarginalR2 = c(0),
                       ConditionalR2 = c(0),
                       effects1 = c("none"),
                       effects2 = c("none"),
                       coefs1 = c(0),
                       coefs2 = c(0),
                       coefsmin1 = c(0),
                       coefsmax1 = c(0),
                       coefsmin2 = c(0),
                       coefsmax2 = c(0),
                       coefsse1 = c(0),
                       coefsse2 = c(0),
)

# fill tibble with stats
for (i in seq_along(best_mod_list)){
  stats_df_veg$MarginalR2[[i]] <- r2(best_mod_list[[i]])[[2]]
  stats_df_veg$ConditionalR2[[i]] <- r2(best_mod_list[[i]])[[1]]
  stats_df_veg$effects1[[i]] <- names(best_mod_list[[i]]$frame)[2]
  stats_df_veg$effects2[[i]] <- names(best_mod_list[[i]]$frame)[3]
  stats_df_veg$coefs1[[i]] <- fixef(best_mod_list[[i]])$cond[2]
  stats_df_veg$coefs2[[i]] <- fixef(best_mod_list[[i]])$cond[3]
  stats_df_veg$coefsmin1[[i]] <- confint(best_mod_list[[i]])[2,1]
  stats_df_veg$coefsmin2[[i]] <- confint(best_mod_list[[i]])[3,1]
  stats_df_veg$coefsmax1[[i]] <- confint(best_mod_list[[i]])[2,2]
  stats_df_veg$coefsmax2[[i]] <- confint(best_mod_list[[i]])[3,2]
  stats_df_veg$coefsse1[[i]] <-
    summary(best_mod_list[[i]])$coefficients$cond[2,2]
  stats_df_veg$coefsse2[[i]] <-
    summary(best_mod_list[[i]])$coefficients$cond[3,2]
}

# drop "Treatment" from fall anthocoridae model (only one fixed effect there)
# use char "NA" to avoid problems with data wrangling below
stats_df_veg2 <- stats_df_veg %>%
  mutate(
    effects2 = case_when(effects2 == "Treatment" ~ "NA",
                         effects2 != "Treatment" ~ effects2),
    coefs2 = case_when(effects2 == "NA" ~ 0,
                         effects2 != "NA" ~ coefs2),
    coefsmin2 = case_when(effects2 == "NA" ~ 0,
                         effects2 != "NA" ~ coefsmin2),
    coefsmax2 = case_when(effects2 == "NA" ~ 0,
                         effects2 != "NA" ~ coefsmax2),
    coefsse2 = case_when(effects2 == "NA" ~ 0,
                         effects2 != "NA" ~ coefsse2)
  )

# Bootstrap ####
# Boostrap data and re-fit models to see whether the top models are robust to
# outliers in the data
## SOURCE ####
source("pred_bootstrap.R", echo = TRUE)

# FIG predator effects ####
# reshape stats table for plotting
figDat <- stats_df %>%
  pivot_longer(effects1:coefsse2,
               names_to = c(".value", "effectRank"),
               names_pattern = "(.*)(.s*)",
               values_to = c("var1, var2, var3, var4, var5, var6")) %>%
  mutate(effects = case_when(effects == "log_AllAph" ~ "logAllAph_NA",
                             effects == "wateringMethod" ~ "wateringMethod_NA",
                             effects == "totalCover" ~ "totalCover_NA",
                             effects == "NA" ~ "NA_NA",
                             # effects == "impermeablePerim" ~ "impermeablePerim_NA",
                             # effects == "naturaAridPerim" ~ "naturalAridPerim_NA",
                             effects != c("log_AllAph",
                                          "wateringMethod",
                                          "totalCover",
                                          "NA") ~ effects)) %>%
  separate(effects, into = c("Effect", "Distance"), sep = "_") %>%
  mutate(Taxon = fct_rev(as_factor(Taxon)),
         Distance = factor(Distance, levels = c("NA",
                                                "sig1",
                                                "sig2",
                                                "sig3",
                                                "sig4",
                                                "sig5",
                                                "const",
                                                "no"))) %>%
  filter(Effect != "NA") %>%
  mutate(moe95 = coefsse*1.96,
         dnum = as.numeric(as.character(recode(Distance,
                                               sig1 = "75",
                                               sig2 = "100",
                                               sig5 = "650",
                                               "NA" = "0",
                                               no = "1000",
                                               const = "999",
                                               sig3 = "350"))),
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
                                totalCover = "Total Cover",
                                divShan = "Land cover diversity",
                                water = "Surface water",
                                # naturalAridPerim = "Desert shrub perimeter",
                                # impermeablePerim = "Impermeable perimeter",
                                # ensure factor order to match palette
                                .ordered = TRUE
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
  ## order factors by L to R appearance
  mutate(Effect = fct_relevel(Effect,
                     "Non-alfalfa agriculture",
                     "Impermeable surfaces",
                     "Weedy cover",
                     "Bare soil",
                     "Alfalfa",
                     "log(Aphid density)",
                     "Land cover diversity",
                     "Flood irrigation"
                     ))

# color palette in correct order
# landcover: Use palette derived from actual colors?
lc_palette_predators <- c(
  "#5b3769", # non-alfalfa  ag
  "grey40", # Impermeable
  "#aecf7a", # weedy cover
  "#801f19", # Bare group
  "#105405", # alfalfa
  "#ff6d00", # log(aphid dens.)
  "#393ebf", # flood irrigation
  "#ff2d55"  # Shannon diversity (land cover)
)

# figDat_boot <- all_boot %>%
#   select(-model) %>%
#   filter(!is.na(effects2)) %>%
#   # mutate_all(~ifelse(is.na(.), "NA", .)) %>%
#   pivot_longer(effects1:coefsse2,
#                names_to = c(".value", "effectRank"),
#                names_pattern = "(.*)(.s*)",
#                values_to = c("var1, var2, var3, var4, var5, var6")) %>%
#   mutate(effects = case_when(effects == "log_AllAph" ~ "logAllAph_NA",
#                              effects == "wateringMethod" ~ "wateringMethod_NA",
#                              effects == "totalCover" ~ "totalCover_NA",
#                              # effects == "impermeablePerim" ~ "impermeablePerim_NA",
#                              # effects == "naturaAridPerim" ~ "naturalAridPerim_NA",
#                              TRUE ~ effects)) %>%
#   separate(effects, into = c("Effect", "Distance"), sep = "_") %>%
#   mutate( Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
#           Taxon = recode_factor(Taxon,
#                                 Anthocoridae = "Anthocor.",
#                                 Arachnida = "Arachnida",
#                                 Coccinellidae = "Coccinell.",
#                                 Geocoris = "Geocoris",
#                                 Ichneumonoidea = "Ichneum."),
#          Distance = factor(Distance, levels = c("NA",
#                                                 "sig1",
#                                                 "sig2",
#                                                 "sig3",
#                                                 "sig4",
#                                                 "sig5",
#                                                 "const",
#                                                 "no"))) %>%
#   mutate(moe95 = coefsse*1.96) %>%
#   # mutate_all(~ifelse(is.na(.), "NA", .)) %>%
#   mutate(dnum = as.numeric(as.character(recode(Distance,
#                                                sig1 = "75",
#                                                sig2 = "100",
#                                                sig5 = "650",
#                                                "NA" = "0",
#                                                no = "1000",
#                                                const = "999",
#                                                sig3 = "350"))),
#          Effect = recode_factor(Effect,
#                                 alfalfa = "Alfalfa",
#                                 natArid = "Desert shrub",
#                                 weedy = "Weedy cover",
#                                 riparian = "Riparian",
#                                 ag = "Non-alfalfa agriculture",
#                                 dirt = "Bare soil",
#                                 impermeable = "Impermeable surfaces",
#                                 logAllAph = "log(Aphid density)",
#                                 wateringMethod = "Flood irrigation",
#                                 totalCover = "Total Cover",
#                                 divShan = "Land cover diversity",
#                                 water = "Surface water",
#                                 # naturalAridPerim = "Desert shrub perimeter",
#                                 # impermeablePerim = "Impermeable perimeter",
#                                 # ensure factor order to match palette
#                                 .ordered = TRUE
#          )) %>%
#   left_join(figDat %>%
#               select(Taxon, Season, Effect, coefsmin, coefsmax) %>%
#               rename("cmin" = coefsmin, "cmax" = coefsmax),
#             by = c("Taxon", "Season", "Effect"))

# make ggplot
ggplot(figDat, aes(y = Taxon,
                   x = coefs,
                   fill = Effect,
                   # color = Effect,
                   moe = moe95)) +
  facet_grid(~Season) +
  stat_confidence_density(height = 0.9,
                          position = position_dodge(0.9),
                          show.legend = TRUE,
                          n = 5000) +
  coord_flip() +
  # scale_y_discrete(drop = FALSE) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(xmin = coefs+coefsse,
                    xmax = coefs-coefsse),
                position = position_dodge(0.9),
                width = 0.2
  ) +
  scale_fill_manual(values = lc_palette_predators) +
  xlim(c(-4,4)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red",
             alpha = 0.5,
             size = 1.2) +
  geom_hline(yintercept = (1:4)+0.5,
             color = "grey70") +
  xlab("Standardized effect size") +
  ylab("Predator Taxon") +
  theme_grey(base_size = 10) +
  theme(#legend.position = c(0.1, 0.1),
    legend.background = element_rect(linetype = 1, color = NA),
    panel.background = element_rect(fill = NA, color = "black"),
    plot.background = element_rect(fill = "white"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 0.8, vjust = 0.95),
    strip.background.x = element_rect(fill = "NA", color = "NA"),
    legend.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin =  margin(r = 0.2, l = -40, t = 0)) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_text(aes(x = 3.7,
                y = Taxon,
                label = cat),
            inherit.aes = F,
            parse = T,
            size = 8 / .pt) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("predator_effects.pdf", width = 18, height = 12, units = "cm", dpi = 600)

# to add data from bootstrapping exercise:
  # geom_jitter(aes(y = Taxon,
  #                x = coefs,
  #                fill = Effect),
  #            data = figDat_boot,
  #            # color = "black",
  #            shape = 21,
  #            position = position_jitterdodge(dodge.width = 0.8),
  #            alpha = 0.5) +
  # geom_text(aes(y = Taxon,
  #               x = coefs,
  #               label = dropped),
  #           data = figDat_boot %>%
  #             filter((coefs < 0 & coefsmax > 0)|(coefs > 0 & coefsmin < 0)),
  #           color = "black",
  #           # position = position_nudge(x = 0.1, y = -0.3),
  #           position = position_dodge(width = 0.8),
  #           size = 2)





## DIVERSION ####
# conclusion: we don't have enough data to explore this.

# look for specific plants associated with Geocoris
# veg_plots %>% View
# look for most variable plants to narrow it down
most_variable <- veg_plots %>%
  group_by(field_id, type) %>%
  summarize(across(SAVE4:UNK_POACEAE5, sum)) %>%
  ungroup() %>%
  group_by(type) %>%
  summarize(across(SAVE4:UNK_POACEAE5, sd)) %>%
  pivot_longer(c(SAVE4:UNK_POACEAE5),
               names_to = "plant",
               values_to = "sd") %>%
  arrange(desc(sd)) %>%
  mutate(plant_type = paste(plant, type, sep = "_"))

ggplot(most_variable %>% slice_head(n = 10),
       aes(fct_rev(fct_reorder(plant_type, sd)), sd)) +
  geom_col()

# Ag: MESA, TRAE, ZEMA, ALCE
# BAHY Margin Bassia hyssopifolia (goosefoot sp. - fivehorn smotherweed)
# SUMO Random Suaeda nigra (chenopodium sp. - bush seepweed)
# ACRE3 Margin Acroptilon repens (thistle sp. - Russian knapweed)
# PLLA Margin Plantago lanceolata
# JUBA Random Juncus balticus (riparian)

geocveg_fa <- veg_plots %>%
  filter(season == "Fall",
         !str_detect(field_id, "Yerington")) %>%
  group_by(field_id, type) %>%
  summarize(across(SAVE4:UNK_POACEAE5, mean), .groups = "keep") %>%
  select(field_id, type, BAHY, SUMO, ACRE3, PLLA, JUBA, ERCI6) %>%
  pivot_wider(names_from = type, values_from = BAHY:ERCI6) %>%
  separate(field_id, into = c("Site", "Field"), sep = " ")

geocData_fa <- df_fa %>%
  left_join(geocveg_fa, by = c("Site", "Field")) %>%
  select(Site, Field, contains("_Margin"), contains("_Random"),
         impermeable_no, wateringMethod, Treatment, Geocoris) %>%
  replace(is.na(.), 0)

geocData_fa %>%
  summarise_all(list(~n_distinct(.))) %>%
  pivot_longer(everything(), names_to = "var")

geocData_fa %>%
  filter(wateringMethod == "Flooding") %>%
  summarise_all(list(~n_distinct(.))) %>%
  pivot_longer(everything(), names_to = "var")

geocData_fa %>%
  filter(wateringMethod == "Sprinklers") %>%
  summarise_all(list(~n_distinct(.))) %>%
  pivot_longer(everything(), names_to = "var")

dotchart(geocData_fa$BAHY_Margin)
dotchart(geocData_fa$SUMO_Margin) # no variation
dotchart(geocData_fa$ACRE3_Margin)
dotchart(geocData_fa$PLLA_Margin) # two levels
dotchart(geocData_fa$JUBA_Margin) # no variation
dotchart(geocData_fa$BAHY_Random) # two levels
dotchart(geocData_fa$SUMO_Random)
dotchart(geocData_fa$ACRE3_Random)
dotchart(geocData_fa$PLLA_Random) # two levels
dotchart(geocData_fa$JUBA_Random) # two levels

geocData_fa %>% arrange(Geocoris) %>% select(Site, Field, Geocoris, ERCI6_Random, ERCI6_Margin)


plant_gmod <- glmmTMB(Geocoris ~
                        BAHY_Margin +
                        SUMO_Margin +
                        ACRE3_Margin +
                        PLLA_Margin +
                        JUBA_Margin +
                        BAHY_Random +
                        SUMO_Random +
                        ACRE3_Random +
                        PLLA_Random +
                        JUBA_Random +
                        # wateringMethod + too few levels!
                        Treatment + (1 | Site/Field),
                      family = "nbinom2",
                      data = geocData_fa,
                      na.action = "na.fail")
clust <- try(makeCluster(getOption("cl.cores", n_cores), type = "PSOCK"))
clusterExport(clust, "geocData_fa")
clusterEvalQ(clust, library(glmmTMB))
plant_gmod_dr <- dredge(plant_gmod, fixed = "cond(Treatment)", m.lim = c(2, 4), trace = 2, cluster = clust)
stopCluster(clust)

# plant_gmod_dr %>% head(10) %>% View

# PLLA_margin?? ACRE3_random???

mod <- get.models(plant_gmod_dr, 1)[[1]]
ggplot(geocData_fa, aes(ACRE3_Random, Geocoris)) +
  geom_point()

plant_gmod <- glmmTMB(Geocoris ~ impermeable_no +
                        wateringMethod + Treatment + (1 | Site/Field),
                      family = "nbinom2",
                      data = geocData_fa)
### end diversion ####

# not robust: fall ara, coc, ich
best_mod_list$best.ara.fa

plot(allEffects(best_mod_list$best.ara.fa))

# choose model to review
review_mod <- best_mod_list$best.geo.fa

# show summary
summary(review_mod) # no random effect variance. essentially equivalent to
# a fixed effects mod. I checked.
# basic effects plots
plot(allEffects(review_mod, residuals = TRUE))

plot(Effect(c("wateringMethod"), review_mod, resid = TRUE))
# try fall mod
reviewFigDf <- ggpredict(review_mod, terms = c("wateringMethod"))
reviewFigDf %>% plot(log.y = F, add.data = T)
# review_mod$frame %>% View

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

source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
influ <- influence_mixed(review_mod)
infIndexPlot(influ)
# View(review_mod$frame)
ggplot(df_fa, aes(weedy_sig3, Arachnida, color = interaction(Site,Field))) +
  geom_point()
ggplot(df_fa, aes(totalCover, Arachnida, color = interaction(Site,Field))) +
  geom_point()
dotchart(df_fa$weedy_sig3, labels = df_fa$Site)
# FIG predator spatial scale ####
# import distance decay dataframe
## TODO incorporate new landcover diversity mods ####
source("distDecay.R", echo = TRUE)

# wrangle figDat to make "nameplates" df
nameplates <- figDat %>%
  # get unique scale value for taxon and season
  mutate(
    combo = paste0(Season, Taxon, dnum),
  ) %>%
  distinct(combo, .keep_all = TRUE) %>%
  # use only landcover vars that have a spatial scale
  filter(dnum > 0) %>%
  # use to manually tune position of geom_label
  mutate(hjust = c(0.2,0.02,-0.08,-0.38,-0.26,
                   0.4,0.08,1.22,1.02,0.1),
         y = c(9,6.5,4,8,1.5,
               8,5.5,3,8,8))  %>% # anth, ara, cocc, geoc, ich
  # make factor for color coding labels
  mutate(col = as.factor(case_when(dnum == 75 ~ "Proximate",
                         dnum == 350 ~ "Intermediate",
                         dnum == 650 ~ "Most distant",
                         dnum == 1000 ~ "None"))
         ) %>%
  # add in missing factor levels
  mutate(col = fct_expand(col, "Near", "Distant", "Constant", "None")) %>%
  # relevel factor order
  mutate(col = fct_relevel(col, "Proximate", "Near", "Intermediate",
                            "Distant", "Most distant", "Constant", "None"))

ggplot(dist_decay_df, aes(distance, value, color = `Decay function`)) +
  geom_line(size = 1) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = dist_palette) +
  scale_fill_manual(values = dist_palette, drop = F, guide = "none") +
  labs(y = "Weighting value") +
  geom_label(data = nameplates, aes(x = dnum,
                             y = y*100,
                             hjust = hjust,
                             label = Taxon,
                             fill = col),
             inherit.aes = FALSE,
             alpha = 0.8,
             size = 8/.pt) +
  facet_wrap(~Season, ncol = 1) +
  theme(strip.background = element_blank(),
        legend.position = "right",
        axis.text = element_text(color = "black"))

ggsave("predator_distance.pdf", width = 18, height = 6, units = "cm", dpi = 600)
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
r2(get.models(tab_nb_allaph_sp_scaled, 1)[[1]]) # good fit 0.85
plot(simulateResiduals(get.models(tab_nb_allaph_sp_scaled, 1)[[1]])) # ok
summary(get.models(tab_nb_allaph_sp_scaled, 1)[[1]])
# tab_nb_allaph_sp_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# basic effects plots
plot(allEffects(get.models(tab_nb_allaph_sp_scaled, 1)[[1]], residuals = TRUE))
# wateringMethod in all top mods, which are close in deltas w/low weights
# -ag_sig1, +alfalfa_sig1, +Cocc

# Acrythosiphon spring
r2(get.models(tab_nb_acy_sp_scaled, 1)[[1]]) # good fit 0.8
plot(simulateResiduals(get.models(tab_nb_acy_sp_scaled, 1)[[1]])) # ok
# tab_nb_acy_sp_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# top mods generally have different landcover factors, but all at _sig1 scale
# watering method still in all top mods
# model average?

# nonacy spring
r2(get.models(tab_nb_nonacy_sp_scaled, 1)[[1]]) # average fit 0.56
plot(simulateResiduals(get.models(tab_nb_nonacy_sp_scaled, 1)[[1]])) # great
# tab_nb_nonacy_sp_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# top mods all include -dirt_no

# Spring Summary NOW OLD
## SEM should include coccinellidae, wateringmethod for sure.
## maybe also include ag_sig1, alfalfa_sig1

### Fall
# AllAph fall
r2(get.models(tab_nb_allaph_fa_scaled, 1)[[1]]) # average fit 0.69
plot(simulateResiduals(get.models(tab_nb_allaph_fa_scaled, 1)[[1]])) # ok
summary(get.models(tab_nb_allaph_fa_scaled, 1)[[1]])
# tab_nb_allaph_fa_scaled %>%
#   tibble %>%
#   slice(1:14) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View

plot(allEffects(get.models(tab_nb_allaph_fa_scaled, 1)[[1]], residuals = TRUE))
# high leverage
dotchart(sort(df_fa$impermeable_sig4), main = "fall impermeable sig4 dotchart")
# df_fa %>% arrange(impermeable_sig4) %>% select(Site, Field, impermeable_sig4, AllAph) %>% View

# Acyrthosiphon fall
r2(get.models(tab_nb_acy_fa_scaled, 1)[[1]]) # pretty good fit 0.7
plot(simulateResiduals(get.models(tab_nb_acy_fa_scaled, 1)[[1]])) # great
# tab_nb_acy_fa_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# +ich in all mods. -geo in top 3 mods. +naturalArid across scales!

# Nonacy fall
r2(get.models(tab_nb_nonacy_fa_scaled, 1)[[1]]) # pretty good fit 0.65
plot(simulateResiduals(get.models(tab_nb_nonacy_fa_scaled, 1)[[1]])) # weird?
# tab_nb_nonacy_fa_scaled %>%
#   tibble %>%
#   slice(1:5) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# +ich in all mods. landcover effects varied in top mods,
# BUT top mod is way better than #2. Includes -impermeable_sig2, +natArid_sig2

# Fall Summary
# strong correlations with ichneumonoidea.
# natural arid benefits aphid abundance. _sig4 best scale for AllAph.

# Review aphid models ####
# spring - #1 mod is inappropriate because it combines wateringmethod and
# water_sig1
sp_best <- get.models(tab_nb_allaph_sp_scaled, 1)[[1]]
fa_best <- get.models(tab_nb_allaph_fa_scaled, 1)[[1]]
summary(sp_best)
summary(fa_best)
# tab_nb_allaph_fa_scaled %>% View
## make aphid modstats table
aph_mods <- list("allaphSP" = sp_best, "allaphFA" = fa_best)
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
                         log_Coccinellidae + shan + rich + totalCover +
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
                         naturalArid_sig4 + log_Ichneumonoidea + shan +
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
# veg_dredge_sp %>%
#   slice(1:15) %>% # can change how inclusive this is
#   select(where(~!all(is.na(.x)))) %>%
#   View
# new best mod!! +richness!!
sp_best_veg <- get.models(veg_dredge_sp, 1)[[1]]
summary(sp_best_veg)
summary(sp_best)
r2(sp_best_veg)
r2(sp_best)
plot(simulateResiduals(sp_best_veg))
plot(fitted(sp_best_veg), residuals(sp_best_veg, type = "pearson"))
plot(df_sp_vd$AllAph, fitted(sp_best_veg))
abline(0, 1)

plot(allEffects(sp_best_veg, resid = TRUE))

## Fall
# veg_dredge_fa %>%
#   slice(1:15) %>% # can change how inclusive this is
#   # select(where(~!all(is.na(.x)))) %>%
#   View
# new best mod!! -shan!!
fa_best_veg <- get.models(veg_dredge_fa, 1)[[1]]
summary(fa_best_veg)
summary(fa_best)
r2(fa_best_veg)
r2(fa_best)
plot(simulateResiduals(fa_best_veg))
plot(fitted(fa_best_veg), residuals(fa_best_veg, type = "pearson"))
plot(df_fa_vd$AllAph, fitted(fa_best_veg))
abline(0, 1)

plot(allEffects(fa_best_veg, resid = TRUE))

# FIG best aphid model effects ####
# check spring effs plot
plot(allEffects(sp_best_veg, resid = TRUE))
plot(Effect(c("rich", "log_Coccinellidae", "wateringMethod"), sp_best_veg, resid = TRUE))
# try spring mod

aphFigDf <- ggpredict(sp_best_veg, terms = c("rich",
                                             "wateringMethod"))
aphFigDfB <- ggpredict(sp_best_veg, terms = c("log_Coccinellidae",
                                             "wateringMethod"))
aphFigDfC <- ggpredict(sp_best_veg, terms = c("log_Coccinellidae",
                                             "rich"))
aphFigDfD <- ggpredict(sp_best_veg, terms = c("log_Coccinellidae",
                                             "rich",
                                             "wateringMethod"))
aphFigDf %>% plot(log.y = T, add.data = T)
aphFigDfB %>% plot(log.y = T, add.data = T)
aphFigDfC %>% plot(log.y = T, add.data = T)
aphFigDfD %>% plot(log.y = T, add.data = T)

# check fall effs plot
plot(allEffects(fa_best_veg, resid = TRUE))
plot(Effect(c("naturalArid_sig4",
              "impermeable_sig4",
              "shan"), fa_best_veg, resid = TRUE))
# try fall mod
aphFigDf2 <- ggpredict(fa_best_veg, terms = c("impermeable_sig4",
                                             "naturalArid_sig4",
                                             "shan"))
# aphFigDf2 %>% plot(log.y = T, add.data = T)

## seems like too many panels. try effects plot similar to predator thing.
## make aphid modstats table
aph_mods_veg <- list("allaphSPVeg" = sp_best_veg, "allaphFAVeg" = fa_best_veg)
# build empty tibble to hold stats
stats_df_aph_veg <- tibble(Taxon = c("AllAph", "AllAph"),
                           Season = c("Spring", "Fall"),
                           MarginalR2 = c(0),
                           ConditionalR2 = c(0),
                           effects1 = c("none"),
                           effects2 = c("none"),
                           effects3 = c("none"),
                           coefs1 = c(0),
                           coefs2 = c(0),
                           coefs3 = c(0),
                           coefsmin1 = c(0),
                           coefsmax1 = c(0),
                           coefsmin2 = c(0),
                           coefsmax2 = c(0),
                           coefsmin3 = c(0),
                           coefsmax3 = c(0),
                           coefsse1 = c(0),
                           coefsse2 = c(0),
                           coefsse3 = c(0)
)

# fill tibble with stats
for (i in seq_along(aph_mods_veg)){
  stats_df_aph_veg$MarginalR2[[i]] <- r2(aph_mods_veg[[i]])[[2]]
  stats_df_aph_veg$ConditionalR2[[i]] <- r2(aph_mods_veg[[i]])[[1]]
  stats_df_aph_veg$effects1[[i]] <- names(aph_mods_veg[[i]]$frame)[2]
  stats_df_aph_veg$effects2[[i]] <- names(aph_mods_veg[[i]]$frame)[3]
  stats_df_aph_veg$effects3[[i]] <- names(aph_mods_veg[[i]]$frame)[4]
  stats_df_aph_veg$coefs1[[i]] <- fixef(aph_mods_veg[[i]])$cond[2]
  stats_df_aph_veg$coefs2[[i]] <- fixef(aph_mods_veg[[i]])$cond[3]
  stats_df_aph_veg$coefs3[[i]] <- fixef(aph_mods_veg[[i]])$cond[4]
  stats_df_aph_veg$coefsmin1[[i]] <- confint(aph_mods_veg[[i]])[2,1]
  stats_df_aph_veg$coefsmin2[[i]] <- confint(aph_mods_veg[[i]])[3,1]
  stats_df_aph_veg$coefsmin3[[i]] <- confint(aph_mods_veg[[i]])[4,1]
  stats_df_aph_veg$coefsmax1[[i]] <- confint(aph_mods_veg[[i]])[2,2]
  stats_df_aph_veg$coefsmax2[[i]] <- confint(aph_mods_veg[[i]])[3,2]
  stats_df_aph_veg$coefsmax3[[i]] <- confint(aph_mods_veg[[i]])[4,2]
  stats_df_aph_veg$coefsse1[[i]] <-
    summary(aph_mods_veg[[i]])$coefficients$cond[2,2]
  stats_df_aph_veg$coefsse2[[i]] <-
    summary(aph_mods_veg[[i]])$coefficients$cond[3,2]
  stats_df_aph_veg$coefsse3[[i]] <-
    summary(aph_mods_veg[[i]])$coefficients$cond[4,2]
}

#### TODO - Export table ####
# plot table for now
stats_df_aph_veg %>%
  tab_df


# FIG aphid effects ####
# reshape stats table for plotting
aphFigDat <- stats_df_aph_veg %>%
  pivot_longer(effects1:coefsse3,
               names_to = c(".value", "effectRank"),
               names_pattern = "(.*)(.s*)",
               values_to = c("var1, var2, var3, var4, var5, var6, var7, var8, var9")) %>%
  mutate(effects = case_when(effects == "log_Coccinellidae" ~ "logCoccinellidae_NA",
                             effects == "log_Ichneumonoidea" ~ "logIchneumonoidea_NA",
                             effects == "wateringMethod" ~ "wateringMethod_NA",
                             effects == "rich" ~ "rich_NA",
                             effects == "shan" ~ "shan_NA",
                             effects == "NA" ~ "NA_NA",
                             effects != c("log_Coccinellidae",
                                          "log_Ichneumonoidea",
                                          "wateringMethod",
                                          "rich",
                                          "shan",
                                          "NA") ~ effects)) %>%
  separate(effects, into = c("Effect", "Distance"), sep = "_") %>%
  mutate(Taxon = fct_rev(as_factor(Taxon)),
         Distance = factor(Distance, levels = c("NA",
                                                "sig1",
                                                "sig2",
                                                "sig3",
                                                "sig4",
                                                "sig5",
                                                "const",
                                                "no"))) %>%
  mutate(moe95 = coefsse*1.96,
         dnum = as.numeric(as.character(recode(Distance,
                                               sig1 = "75",
                                               sig4 = "650",
                                               "NA" = "0",
                                               no = "1000",
                                               sig3 = "350"))),
         Effect = recode_factor(Effect,
                                naturalArid = "Desert shrub,\n\"Distant\" range",
                                impermeable = "Impermeable surfaces,\n\"Distant\" range",
                                logCoccinellidae = "log(Coccinellidae density)",
                                logIchneumonoidea = "log(Ichneumonoidea density)",
                                wateringMethod = "Flood irrigation",
                                shan = "Shannon Diversity\n(in field margins)",
                                rich = "Plant species richness\n(in field margins)",
                                # ensure factor order to match palette
                                .ordered = TRUE
         ),
         Season = factor(Season, levels = c(Spring = "Spring", Fall = "Fall")),
         Taxon = recode_factor(Taxon,
                               AllAph = "Aphid density"),
         roundM = round(MarginalR2, 3),
         roundC = round(ConditionalR2, 3),
         quo = "{R^2}[M]~'='",
         quoC = "{R^2}[C]~'='",
         cat = paste0(quo, "~'", roundM, "'"),
         catC = paste0(quoC, "~'", roundC, "'"),
         expr2 = list(bquote("{R^2}[M]~'='"~.(roundM))),
         expr2C = list(bquote("{R^2}[C]~'='"~.(roundC))),
         expr = paste0(round(MarginalR2,2)),
         exprC = paste0(round(ConditionalR2,2)))

ggplot(aphFigDat, aes(y = Effect,
                   x = coefs,
                   # fill = Effect,
                   # color = Effect,
                   moe = moe95)) +
  facet_wrap(~Season, ncol = 1, scales = "free_y") +
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
  # scale_fill_manual(values = lc_palette_experimental) +
  xlim(c(-2,2)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red",
             alpha = 0.5) +
  xlab("Standardized effect size") +
  ylab(NULL) +
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
    axis.text = element_text(color = "black", size = 6),
    axis.title.x = element_text(size = 6),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box.margin =  margin(r = 0.2, l = -40, t = 0)) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_text(aes(x = -1.5,
                y = 1.5,
                label = cat),
            inherit.aes = F,
            parse = T,
            size = 6/.pt) +
  geom_text(aes(x = -1.5,
                y = 0.8,
                label = catC),
            inherit.aes = F,
            parse = T,
            size = 6/.pt)

ggsave("aphid_effects.pdf", width = 8.5, height = 6, units = "cm", dpi = 600)


df_fa_vd %>%
  ggplot(aes(Site, impermeable_sig4, color = Field)) +
  geom_jitter()
df_fa %>%
  ggplot(aes(alfalfa_sig4, impermeable_sig4, color = Site)) +
  geom_jitter()+
  geom_smooth(method = "lm")
df_fa %>%
  ggplot(aes(log_Ichneumonoidea, impermeable_sig4, color = Site)) +
  geom_jitter()

df_fa_vd %>%
  ggplot(aes(impermeable_sig4, log_AllAph, color = shan)) +
  geom_point()

fa_best_veg %>% summary()

dotplot(df_fa$impermeable_sig5)
dotplot(df_fa$impermeable_sig4)
dotplot(df_fa$impermeable_sig4)
dotplot(df_fa$impermeable_sig3)
dotplot(df_fa$impermeable_sig2)
dotplot(df_fa$impermeable_sig1)

## FIG watering method ####

library(ggbeeswarm)


flood_mean <- pull(df_sp %>% filter(wateringMethod == "Flooding"),
     log_AllAph) %>% mean

sprinklers_mean <- pull(df_sp %>% filter(wateringMethod == "Sprinklers"),
     log_AllAph) %>% mean

ggblue <- scales::hue_pal()(2)[2]
ggred <- scales::hue_pal()(2)[1]

p <-ggplot(df_sp, aes(wateringMethod, log_AllAph, fill = wateringMethod))+
  geom_boxplot() +
  geom_quasirandom(bandwidth = 0.3) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA, color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  labs(y = "log(Aphid density)",
       x = "Watering method")

q <- ggplot(df_sp, aes(water_no, log_AllAph, color = wateringMethod))+
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  # geom_hline(yintercept = flood_mean,
  #            color = ggblue) +
  # geom_hline(yintercept = sprinklers_mean,
  #            color = ggred) +
  # geom_text(aes(-0.9, flood_mean+0.1),
  #           label = "Flooding - mean",
  #           color = ggblue,
  #           size = 8/.pt) +
  # geom_text(aes(1, sprinklers_mean+0.1),
  #           label = "Sprinklers - mean",
  #           color = ggred,
  #           size = 8/.pt) +
  theme_grey(base_size = 10) +
  theme(legend.position = c(0.9, 0.15),
        legend.background = element_rect(linetype = 1, color = NA),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        # axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.box.margin =  margin(r = 0.2, l = -40, t = 0),
        legend.key = element_blank()) +
  labs(x = "Surface water, no distance weighting",
       y = "log(Aphid density)")


p
q


ggsave(file = "flood_effect.pdf",
       plot = marrangeGrob(list(q, p), ncol = 2, nrow = 1, top = NULL),
       width = 18, height = 8.5, units = "cm", dpi = 600)



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
  left_join(diff_data_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-") %>%
         # diffCoccinellidae > 0) %>%
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) %>%
  select(logAllAph, logCoccinellidae, Coccinellidae, # dump other cols
         diffCoccinellidae,
         wateringMethod_Flooding,
         ag_sig1, weedy_sig1, dirt_sig1,
         AllAph, Treatment, rich,
         Treatment_Sham, Site, Field) %>%
  mutate(diffCoccinellidae =
           case_when(Treatment_Sham == 1 ~ diffCoccinellidae,
                     Treatment_Sham == 0 ~ 0)) %>%
  mutate(across(.cols = -c(Treatment, Site, Field, Treatment_Sham,
                           wateringMethod_Flooding, diffCoccinellidae),
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
  left_join(diff_data_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-",
         diffCoccinellidae > 0) %>% # filter here!
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) %>%
  select(logAllAph, logCoccinellidae, Coccinellidae, # dump other cols
         diffCoccinellidae, rich,
         wateringMethod_Flooding,
         ag_sig1, weedy_sig1, dirt_sig1,
         AllAph, Treatment,
         Treatment_Sham, Site, Field) %>%
  mutate(diffCoccinellidae =
           case_when(Treatment_Sham == 1 ~ diffCoccinellidae,
                     Treatment_Sham == 0 ~ 0)) %>%
  mutate(across(.cols = -c(Treatment, Site, Field,
                           wateringMethod_Flooding, rich,
                           Treatment_Sham, diffCoccinellidae),
                .fns = ~(./sd(.)))) # scale, don't center
## specify model (with no constraints)
mod_specA <- "
  # direct effects
    logCoccinellidae ~ w*weedy_sig1 + d*dirt_sig1

    logAllAph ~ f*wateringMethod_Flooding + a*rich + t*Treatment_Sham

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
# Trtmnt_Shm (t) Est: -0.847    Std.err: 0.126   z: -6.732    P: 0.000

lavaanPlot(model = mod_fitA, coefs = TRUE, covs = F, stand = FALSE)


## 3b. refit to full dataset with Trt_Sham coef "fixed"
mod_specB <- "
  # direct effects
    logCoccinellidae ~ w*weedy_sig1 + d*dirt_sig1

    logAllAph ~ f*wateringMethod_Flooding + a*rich + t*Treatment_Sham

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
   t == -0.847
"

# ML is default estimator
mod_fit1B <- sem(mod_specB, data = lavaan_df1)
summary(mod_fit1B)

lavaanPlot(model = mod_fit1B, coefs = TRUE, covs = FALSE, stand = FALSE)
lavaanPlot(model = mod_fit1B, coefs = TRUE, covs = TRUE, stand = FALSE)


lavResiduals(mod_fit1B)

### fall lavaan ####
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
  left_join(diff_data_wide) %>%
  filter(Season == "Fall",
         Treatment != "Pre-") %>%
  # diffCoccinellidae > 0) %>%
  mutate(logAllAph = log(AllAph + 1),
         logCoccinellidae = log(Coccinellidae + 1)) %>%
  select(logAllAph, impermeable_sig4, naturalArid_sig4, shan,
         AllAph, Treatment,
         Treatment_Sham, Site, Field) %>%
  mutate(across(.cols = -c(Treatment, Site, Field, Treatment_Sham,
                           shan),
                .fns = ~as.vector(scale(.x))))
# alternative: ~(. -mean(.)/ sd(.))
## specify model (with no constraints)
mod_specA <- "
  # direct effects

    logAllAph ~ i*impermeable_sig4 + n*naturalArid_sig4 + t*Treatment_Sham

  # mediator
    # ?

  # indirect effects
    # ?
    # multiply coefs here - use :=

  # total effects
    # add coefs here - use :=

  # covariance
    # add covs here

  # constraints
   # none yet!
"

mod_fitA <- sem(mod_specA, data = lavaan_df1)

summary(mod_fitA)
## WHAT IS THIS BELOW? ##########
# cant figure out where this came from
# Trt_Sham effect is:
# Trtmnt_Shm (t) Est: -0.758    Std.err: 0.133   z: -5.708    P: 0.000
#################################

lavaanPlot(model = mod_fitA, coefs = TRUE, covs = F, stand = FALSE)
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

# what is the deal with the diffCoccinellidae?
