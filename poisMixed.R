## this script must be sourced by data_analysis ##

# glmmTMB with nested Site/Field random effects (where possible)
# nbinom2 distribution
# scaled and centered explanatory variables

# not following style conventions to save space.
# this code is necessarily repetitive.

# set global NA action
options(na.action = 'na.fail')
# make it parallel
clust <- try(makeCluster(getOption("cl.cores", n.cores), type = 'PSOCK'))
clusterExport(clust, 'dfSp')
clusterExport(clust, 'dfFa')
clusterEvalQ(clust, library(glmmTMB))

# Anthocoridae ####
## Spring ####
gMod.sig1 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig2 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig3 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig4 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig5 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.const <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = poisson())
gMod.no <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.anth.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site), # nested random effects NOT fitted - convergence error. Simplified to site only effect.
                     data = dfFa, family = poisson())
gMod.sig2 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig3 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig4 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig5 <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.const <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = poisson()) # convergence warning
gMod.no <- glmmTMB(Anthocoridae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.anth.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)



# Arachnida ####
## Spring ####
gMod.sig1 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig2 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig3 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson()) # convergence warning
gMod.sig4 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig5 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.const <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = poisson())
gMod.no <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.ara.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # Random effects NOT fitted (convergence error)
                     data = dfFa, family = poisson())
gMod.sig2 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig3 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig4 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig5 <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.const <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = poisson())
gMod.no <- glmmTMB(Arachnida ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.ara.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)


# Coccinellidae ####
## Spring ####
gMod.sig1 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig2 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig3 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig4 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig5 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.const <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = poisson())
gMod.no <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.cocc.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site), # nested random effects NOT fitted (convergence error). Simplified to include site only.
                     data = dfFa, family = poisson())
gMod.sig2 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig3 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig4 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig5 <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.const <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = poisson())
gMod.no <- glmmTMB(Coccinellidae ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.cocc.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)





# Geocoris ####
## Spring ####
gMod.sig1 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig2 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig3 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig4 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig5 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.const <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = poisson())
gMod.no <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.geo.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig2 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig3 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig4 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig5 <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.const <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = poisson()) # convergence warning
gMod.no <- glmmTMB(Geocoris ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = poisson()) # convergence warning

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.geo.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)






# Ichneumonoidea ####
## Spring ####
gMod.sig1 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig2 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig3 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig4 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.sig5 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = poisson())
gMod.const <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = poisson())
gMod.no <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.ich.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # Random effects NOT fitted (convergence error)
                     data = dfFa, family = poisson())
gMod.sig2 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig3 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson()) # convergence warning
gMod.sig4 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.sig5 <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = poisson())
gMod.const <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = poisson())
gMod.no <- glmmTMB(Ichneumonoidea ~ Treatment + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = poisson())

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.pois.ich.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

stopCluster(clust)

# put all tabs in lists
pois.scaled <- list(
  'tab.pois.anth.sp.scaled'=tab.pois.anth.sp.scaled,
  'tab.pois.anth.fa.scaled'=tab.pois.anth.fa.scaled,
  'tab.pois.ara.sp.scaled'=tab.pois.ara.sp.scaled,
  'tab.pois.ara.fa.scaled'=tab.pois.ara.fa.scaled,
  'tab.pois.cocc.sp.scaled'=tab.pois.cocc.sp.scaled,
  'tab.pois.cocc.fa.scaled'=tab.pois.cocc.fa.scaled,
  'tab.pois.geo.sp.scaled'=tab.pois.geo.sp.scaled,
  'tab.pois.geo.fa.scaled'=tab.pois.geo.fa.scaled,
  'tab.pois.ich.sp.scaled'=tab.pois.ich.sp.scaled,
  'tab.pois.ich.fa.scaled'=tab.pois.ich.fa.scaled
)

# clean env
rm(
  tab.pois.anth.sp.scaled,
  tab.pois.anth.fa.scaled,
  tab.pois.ara.sp.scaled,
  tab.pois.ara.fa.scaled,
  tab.pois.cocc.sp.scaled,
  tab.pois.cocc.fa.scaled,
  tab.pois.geo.sp.scaled,
  tab.pois.geo.fa.scaled,
  tab.pois.ich.sp.scaled,
  tab.pois.ich.fa.scaled,
  sig1.dredge,
  sig2.dredge,
  sig3.dredge,
  sig4.dredge,
  sig5.dredge,
  const.dredge,
  no.dredge,
  gMod.sig1,
  gMod.sig2,
  gMod.sig3,
  gMod.sig4,
  gMod.sig5,
  gMod.const,
  gMod.no
)
