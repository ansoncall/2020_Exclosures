# nb mixed aph

## this script must be sourced by data_analysis ##

# glmmTMB with nested Site / Field random effects (where possible)
# nbinom2 distribution
# scaled and centered explanatory variables

# not following style conventions to save space_
# this code is necessarily repetitive_

# set global NA action
options(na.action = "na.fail")
# make it parallel
clust <- try(makeCluster(getOption("cl.cores", n_cores), type = "PSOCK"))
clusterExport(clust, "df_sp")
clusterExport(clust, "df_fa")
clusterEvalQ(clust, library(glmmTMB))

# AllAph ####
## Spring ####
gmod_sig1 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig2 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig3 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig4 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig5 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_const <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = "nbinom2")
gmod_no <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_allaph_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Geocoris + log_Ichneumonoidea, # nested random effects not fitted
                     data = df_fa, family = "nbinom2")
gmod_sig2 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig3 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig4 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig5 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_const <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = "nbinom2")
gmod_no <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_allaph_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

# Acyrthosiphon ####
## Spring ####
gmod_sig1 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig2 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig3 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig4 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig5 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_const <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = "nbinom2")
gmod_no <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_acy_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Geocoris + log_Ichneumonoidea, # nested random effects not fitted
                     data = df_fa, family = "nbinom2")
gmod_sig2 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig3 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig4 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig5 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_const <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = "nbinom2")
gmod_no <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_acy_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

# NonAcy ####
## Spring ####
gmod_sig1 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig2 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig3 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig4 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_sig5 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = "nbinom2")
gmod_const <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = "nbinom2")
gmod_no <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_nonacy_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 + divShan_sig1 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea, # nested random effects not fitted
                     data = df_fa, family = "nbinom2")
gmod_sig2 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 + divShan_sig2 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig3 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 + divShan_sig3 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig4 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 + divShan_sig4 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_sig5 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 + divShan_sig5 + # landcover effects
                       + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = "nbinom2")
gmod_const <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + div_const + divShan_const + # landcover effects
                        + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = "nbinom2")
gmod_no <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + div_no + divShan_no + # landcover effects
                     + log_Anthocoridae + log_Arachnida + log_Coccinellidae + log_Geocoris + log_Ichneumonoidea + # predator effects
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = "nbinom2")

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 4), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_nb_nonacy_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)
