# Compare vegetation - Predators

# This script takes the best model from the model selection table, and creates
# the same model for a reduced dataset (-Yerington) with additional parameters
# from the vegetation survey data. It then dredges this new global model and
# returns the best model from the new dredge.

# Spring ####
## Anthocoridae
best_mod_list$best.ant.sp
global.ant.sp.vd <- glmmTMB(Anthocoridae ~ Treatment+ag_sig1 +
                              impermeable_sig1 + (1|Site:Field) + shan + rich +
                              totalCover, data = df_sp_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ant.sp.vd.tab <- dredge(global.ant.sp.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ant.sp.vd <- get.models(ant.sp.vd.tab, 1)[[1]]

## Arachnida
best_mod_list$best.ara.sp
global.ara.sp.vd <- glmmTMB(Arachnida ~ Treatment + dirt_sig1 + weedy_sig1 +
                              (1|Site:Field) + shan + rich + totalCover,
                            data = df_sp_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ara.sp.vd.tab <- dredge(global.ara.sp.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ara.sp.vd <- get.models(ara.sp.vd.tab, 1)[[1]]

## Coccinellidae
best_mod_list$best.coc.sp
global.coc.sp.vd <- glmmTMB(Coccinellidae ~ Treatment + dirt_sig1 + weedy_sig1 +
                              (1|Site:Field) + shan + rich + totalCover,
                            data = df_sp_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
coc.sp.vd.tab <- dredge(global.coc.sp.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.coc.sp.vd <- get.models(coc.sp.vd.tab, 1)[[1]]

## Ichneumonoidea
best_mod_list$best.ich.sp
global.ich.sp.vd <- glmmTMB(Ichneumonoidea ~ Treatment + impermeable_sig1 +
                              log(AllAph + 1) + (1|Site:Field) + shan + rich +
                              totalCover,
                            data = df_sp_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ich.sp.vd.tab <- dredge(global.ich.sp.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ich.sp.vd <- get.models(ich.sp.vd.tab, 1)[[1]]

# Fall ####
# Anthocoridae
best_mod_list$best.ant.fa
global.ant.fa.vd <- glmmTMB(Anthocoridae ~ Treatment + alfalfa_sig1 +
                              dirt_sig1 + (1|Site:Field) + shan + rich +
                              totalCover,
                            data = df_fa_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ant.fa.vd.tab <- dredge(global.ant.fa.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ant.fa.vd <- get.models(ant.fa.vd.tab, 1)[[1]]

# Arachnida
best_mod_list$best.ara.fa
global.ara.fa.vd <- glmmTMB(Arachnida ~ Treatment + weedy_sig3 + (1|Site:Field)
                            + shan + rich + totalCover,
                            data = df_fa_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ara.fa.vd.tab <- dredge(global.ara.fa.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ara.fa.vd <- get.models(ara.fa.vd.tab, 1)[[1]]
# must drop wateringMethod here because all sites are flooded.

# Coccinellidae
best_mod_list$best.coc.fa
global.coc.fa.vd <- glmmTMB(Coccinellidae ~ Treatment + ag_no + alfalfa_no +
                              (1|Site:Field) + shan + rich + totalCover,
                            data = df_fa_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
coc.fa.vd.tab <- dredge(global.coc.fa.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.coc.fa.vd <- get.models(coc.fa.vd.tab, 1)[[1]]

# Ichneumonoidea
best_mod_list$best.ich.fa
global.ich.fa.vd <- glmmTMB(Ichneumonoidea ~ Treatment + ag_sig1 +
                              log(AllAph + 1) + (1|Site:Field) + shan + rich +
                              totalCover,
                            data = df_fa_vd,
                            family = 'nbinom2',
                            na.action = "na.fail")
ich.fa.vd.tab <- dredge(global.ich.fa.vd,
                        m.lim = c(0,3),
                        fixed = 'cond(Treatment)',
                        trace = 2)
best.ich.fa.vd <- get.models(ich.fa.vd.tab, 1)[[1]]
