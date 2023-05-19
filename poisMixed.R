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

# Anthocoridae ####
## Spring ####
gmod_sig1 <- glmmTMB(Anthocoridae ~
                       Treatment + log_AllAph + wateringMethod +

                       # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig2 <- glmmTMB(Anthocoridae ~
                       Treatment + log_AllAph + wateringMethod +
                       # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig3 <- glmmTMB(Anthocoridae ~
                       Treatment + log_AllAph + wateringMethod +
                       # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig4 <- glmmTMB(Anthocoridae ~
                       Treatment + log_AllAph + wateringMethod +
                       # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig5 <- glmmTMB(Anthocoridae ~
                       Treatment + log_AllAph + wateringMethod +
                       # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       # landcover effects
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_const <- glmmTMB(Anthocoridae ~
                        Treatment + log_AllAph + wateringMethod +
                        # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        # landcover effects
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = poisson())
gmod_no <- glmmTMB(Anthocoridae ~
                     Treatment + log_AllAph + wateringMethod +
                     # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +


                     # landcover effects
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_anth_sp_scaled <- rbind(sig1_dredge,
                                 sig2_dredge,
                                 sig3_dredge,
                                 sig4_dredge,
                                 sig5_dredge,
                                 const_dredge,
                                 no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Anthocoridae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1, + div_sig2 +
                       divShan_sig2,
                     # nested random effects NOT fitted
                     data = df_fa, family = poisson())
gmod_sig2 <- glmmTMB(Anthocoridae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig3 <- glmmTMB(Anthocoridae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig4 <- glmmTMB(Anthocoridae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig5 <- glmmTMB(Anthocoridae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_const <- glmmTMB(Anthocoridae ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = poisson()) # convergence warning
gmod_no <- glmmTMB(Anthocoridae ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_anth_fa_scaled <- rbind(sig1_dredge,
                                 sig2_dredge,
                                 sig3_dredge,
                                 sig4_dredge,
                                 sig5_dredge,
                                 const_dredge,
                                 no_dredge)

# Arachnida ####
## Spring ####
gmod_sig1 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig2 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig3 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig5 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_const <- glmmTMB(Arachnida ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const,
                      # (1 | Site / Field), # nested random effects not fitted
                      data = df_sp, family = poisson())
gmod_no <- glmmTMB(Arachnida ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_ara_sp_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1,
                     # Random effects NOT fitted (convergence error)
                     data = df_fa, family = poisson())
gmod_sig2 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig4 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig5 <- glmmTMB(Arachnida ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_const <- glmmTMB(Arachnida ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = poisson())
gmod_no <- glmmTMB(Arachnida ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_ara_fa_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

# Coccinellidae ####
## Spring ####
gmod_sig1 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig2 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig3 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig4 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig5 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_const <- glmmTMB(Coccinellidae ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = poisson())
gmod_no <- glmmTMB(Coccinellidae ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_cocc_sp_scaled <- rbind(sig1_dredge,
                                 sig2_dredge,
                                 sig3_dredge,
                                 sig4_dredge,
                                 sig5_dredge,
                                 const_dredge,
                                 no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1,
                     # nested random effects NOT fitted (convergence error)_
                     data = df_fa, family = poisson())
gmod_sig2 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig5 <- glmmTMB(Coccinellidae ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_const <- glmmTMB(Coccinellidae ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = poisson())
gmod_no <- glmmTMB(Coccinellidae ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_cocc_fa_scaled <- rbind(sig1_dredge,
                                 sig2_dredge,
                                 sig3_dredge,
                                 sig4_dredge,
                                 sig5_dredge,
                                 const_dredge,
                                 no_dredge)

# Geocoris ####
## Spring ####
gmod_sig1 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig2 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig3 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3,
                     # (1 | Site / Field), # nested random effects not fitted
                     data = df_sp, family = poisson())
gmod_sig4 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig5 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_const <- glmmTMB(Geocoris ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = poisson())
gmod_no <- glmmTMB(Geocoris ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no,
                   # (1 | Site / Field), # nested random effects not fitted
                   data = df_sp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_geo_sp_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig2 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig4 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig5 <- glmmTMB(Geocoris ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_const <- glmmTMB(Geocoris ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = poisson()) # convergence warning
gmod_no <- glmmTMB(Geocoris ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = poisson()) # convergence warning

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_geo_fa_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

# Ichneumonoidea ####
## Spring ####
gmod_sig1 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1 + div_sig1 +
                       divShan_sig1 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig2 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig3 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig4 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_sig5 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_sp, family = poisson())
gmod_const <- glmmTMB(Ichneumonoidea ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_sp, family = poisson())
gmod_no <- glmmTMB(Ichneumonoidea ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_sp, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_ich_sp_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

## Fall ####
gmod_sig1 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 + ag_sig1 +
                       impermeable_sig1 + weedy_sig1 + water_sig1,
                     # Random effects NOT fitted (convergence error)
                     data = df_fa, family = poisson())
gmod_sig2 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 + ag_sig2 +
                       impermeable_sig2 + weedy_sig2 + water_sig2 + div_sig2 +
                       divShan_sig2 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig3 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 + ag_sig3 +
                       impermeable_sig3 + weedy_sig3 + water_sig3 + div_sig3 +
                       divShan_sig3 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson()) # convergence warning
gmod_sig4 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 + ag_sig4 +
                       impermeable_sig4 + weedy_sig4 + water_sig4 + div_sig4 +
                       divShan_sig4 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_sig5 <- glmmTMB(Ichneumonoidea ~
                       # non-landcover effects
                       Treatment + log_AllAph + wateringMethod +
                       # landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 + ag_sig5 +
                       impermeable_sig5 + weedy_sig5 + water_sig5 + div_sig5 +
                       divShan_sig5 +
                       (1 | Site / Field), # nested random effects
                     data = df_fa, family = poisson())
gmod_const <- glmmTMB(Ichneumonoidea ~
                        # non-landcover effects
                        Treatment + log_AllAph + wateringMethod +
                        # landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +
                        ag_const + impermeable_const + weedy_const +
                        water_const + div_const + divShan_const +
                        (1 | Site / Field), # nested random effects
                      data = df_fa, family = poisson())
gmod_no <- glmmTMB(Ichneumonoidea ~
                     # non-landcover effects
                     Treatment + log_AllAph + wateringMethod +
                     # landcover effects
                     alfalfa_no + naturalArid_no + dirt_no + ag_no +
                     impermeable_no + weedy_no + water_no + div_no +
                     divShan_no +
                     (1 | Site / Field), # nested random effects
                   data = df_fa, family = poisson())

# dredging
sig1_dredge <- dredge(gmod_sig1,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig2_dredge <- dredge(gmod_sig2,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig3_dredge <- dredge(gmod_sig3,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig4_dredge <- dredge(gmod_sig4,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
sig5_dredge <- dredge(gmod_sig5,
                      m.lim = c(0, 3),
                      fixed = "cond(Treatment)",
                      trace = 2,
                      cluster = clust)
const_dredge <- dredge(gmod_sigconst,
                       m.lim = c(0, 3),
                       fixed = "cond(Treatment)",
                       trace = 2,
                       cluster = clust)
no_dredge <- dredge(gmod_signo,
                    m.lim = c(0, 3),
                    fixed = "cond(Treatment)",
                    trace = 2,
                    cluster = clust)

# rbind and recalc AIC
tab_pois_ich_fa_scaled <- rbind(sig1_dredge,
                                sig2_dredge,
                                sig3_dredge,
                                sig4_dredge,
                                sig5_dredge,
                                const_dredge,
                                no_dredge)

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
  gmod_sigconst,
  gmod_no
)
