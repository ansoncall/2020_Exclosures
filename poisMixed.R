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
clusterExport(clust, "dfSp")
clusterExport(clust, "dfFa")
clusterEvalQ(clust, library(glmmTMB))

# Anthocoridae ####
## Spring ####
gmod_sig1 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig2 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig3 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig4 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig5 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_const <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfSp, family = poisson())
gmod_no <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_anth_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects NOT fitted
                     data = dfFa, family = poisson())
gmod_sig2 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig3 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig4 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig5 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_const <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfFa, family = poisson()) # convergence warning
gmod_no <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_anth_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)



# Arachnida ####
## Spring ####
gmod_sig1 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig2 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig3 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig5 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_const <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfSp, family = poisson())
gmod_no <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_ara_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # Random effects NOT fitted (convergence error)
                     data = dfFa, family = poisson())
gmod_sig2 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig4 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig5 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_const <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfFa, family = poisson())
gmod_no <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_ara_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)


# Coccinellidae ####
## Spring ####
gmod_sig1 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig2 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig3 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig4 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig5 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_const <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfSp, family = poisson())
gmod_no <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_cocc_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects NOT fitted (convergence error)_
                     data = dfFa, family = poisson())
gmod_sig2 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig5 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_const <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfFa, family = poisson())
gmod_no <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_cocc_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)





# Geocoris ####
## Spring ####
gmod_sig1 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig2 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig3 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig4 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig5 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_const <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfSp, family = poisson())
gmod_no <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_geo_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig2 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig4 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig5 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_const <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfFa, family = poisson()) # convergence warning
gmod_no <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfFa, family = poisson()) # convergence warning

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_geo_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)






# Ichneumonoidea ####
## Spring ####
gmod_sig1 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig2 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig3 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig4 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_sig5 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfSp, family = poisson())
gmod_const <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfSp, family = poisson())
gmod_no <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_ich_sp_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # Random effects NOT fitted (convergence error)
                     data = dfFa, family = poisson())
gmod_sig2 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_sig5 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = dfFa, family = poisson())
gmod_const <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const + ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = dfFa, family = poisson())
gmod_no <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig2_dredge <- dredge(gmod_sig2, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig3_dredge <- dredge(gmod_sig3, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig4_dredge <- dredge(gmod_sig4, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
sig5_dredge <- dredge(gmod_sig5, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
const_dredge <- dredge(gmod_const, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)
no_dredge <- dredge(gmod_no, m.lim = c(0, 3), fixed = "cond(Treatment)", trace = 2, cluster = clust)

# rbind and recalc AIC
tab_pois_ich_fa_scaled <- rbind(sig1_dredge, sig2_dredge, sig3_dredge, sig4_dredge, sig5_dredge, const_dredge, no_dredge)

stopCluster(clust)

# put all tabs in lists
pois_scaled <- list(
  "tab_pois_anth_sp_scaled" = tab_pois_anth_sp_scaled,
  "tab_pois_anth_fa_scaled" = tab_pois_anth_fa_scaled,
  "tab_pois_ara_sp_scaled" = tab_pois_ara_sp_scaled,
  "tab_pois_ara_fa_scaled" = tab_pois_ara_fa_scaled,
  "tab_pois_cocc_sp_scaled" = tab_pois_cocc_sp_scaled,
  "tab_pois_cocc_fa_scaled" = tab_pois_cocc_fa_scaled,
  "tab_pois_geo_sp_scaled" = tab_pois_geo_sp_scaled,
  "tab_pois_geo_fa_scaled" = tab_pois_geo_fa_scaled,
  "tab_pois_ich_sp_scaled" = tab_pois_ich_sp_scaled,
  "tab_pois_ich_fa_scaled" = tab_pois_ich_fa_scaled
)

# clean env
rm(
  tab_pois_anth_sp_scaled,
  tab_pois_anth_fa_scaled,
  tab_pois_ara_sp_scaled,
  tab_pois_ara_fa_scaled,
  tab_pois_cocc_sp_scaled,
  tab_pois_cocc_fa_scaled,
  tab_pois_geo_sp_scaled,
  tab_pois_geo_fa_scaled,
  tab_pois_ich_sp_scaled,
  tab_pois_ich_fa_scaled,
  sig1_dredge,
  sig2_dredge,
  sig3_dredge,
  sig4_dredge,
  sig5_dredge,
  const_dredge,
  no_dredge,
  gmod_sig1,
  gmod_sig2,
  gmod_sig3,
  gmod_sig4,
  gmod_sig5,
  gmod_const,
  gmod_no
)
