## this script must be sourced by data_analysis ##

# glmmTMB with nested Site/Field random effects (where possible)
# nbinom2 distribution
# scaled and centered explanatory variables

# not following style conventions to save space.
# this code is necessarily repetitive.

# set global NA action
options(na.action = 'na.fail')
# make it parallel
clust <- try(makeCluster(getOption("cl.cores", 12), type = 'PSOCK'))
clusterExport(clust, 'dfSp')
clusterExport(clust, 'dfFa')
clusterEvalQ(clust, library(glmmTMB))

# Anthocoridae ####
## Spring ####
gMod.sig1 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.anth.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects not fitting.
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Anthocoridae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.anth.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)



# Arachnida ####
## Spring ####
gMod.sig1 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.ara.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects not fitted.
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Arachnida ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.ara.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)


# Coccinellidae ####
## Spring ####
gMod.sig1 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.cocc.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects not fitted
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Coccinellidae ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.cocc.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)





# Geocoris ####
## Spring ####
gMod.sig1 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.geo.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Geocoris ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.geo.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)






# Ichneumonoidea ####
## Spring ####
gMod.sig1 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.ich.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1, # landcover effects
                       # nested random effects not fitted.
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Ichneumonoidea ~ offset(Area) + log_AllAph + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 3), fixed = 'cond(offset(Area))', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.ich.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

stopCluster(clust)

# put all tabs in lists
nb.scaled <- list(
  'tab.nb.anth.sp.scaled'=tab.nb.anth.sp.scaled,
  'tab.nb.anth.fa.scaled'=tab.nb.anth.fa.scaled,
  'tab.nb.ara.sp.scaled'=tab.nb.ara.sp.scaled,
  'tab.nb.ara.fa.scaled'=tab.nb.ara.fa.scaled,
  'tab.nb.cocc.sp.scaled'=tab.nb.cocc.sp.scaled,
  'tab.nb.cocc.fa.scaled'=tab.nb.cocc.fa.scaled,
  'tab.nb.geo.sp.scaled'=tab.nb.geo.sp.scaled,
  'tab.nb.geo.fa.scaled'=tab.nb.geo.fa.scaled,
  'tab.nb.ich.sp.scaled'=tab.nb.ich.sp.scaled,
  'tab.nb.ich.fa.scaled'=tab.nb.ich.fa.scaled
)

# clean env
rm(
  tab.nb.anth.sp.scaled,
  tab.nb.anth.fa.scaled,
  tab.nb.ara.sp.scaled,
  tab.nb.ara.fa.scaled,
  tab.nb.cocc.sp.scaled,
  tab.nb.cocc.fa.scaled,
  tab.nb.geo.sp.scaled,
  tab.nb.geo.fa.scaled,
  tab.nb.ich.sp.scaled,
  tab.nb.ich.fa.scaled,
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
