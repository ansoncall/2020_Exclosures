# nb mixed aph

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

# AllAph ####
## Spring ####
gMod.sig1 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.allaph.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Geocoris+1) + log(Ichneumonoidea+1), # nested random effects not fitted
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(AllAph ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(AllAph ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.allaph.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

# Acyrthosiphon ####
## Spring ####
gMod.sig1 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.acy.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Geocoris+1) + log(Ichneumonoidea+1), # nested random effects not fitted
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(Acyrthosiphon ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(Acyrthosiphon ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.acy.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

# NonAcy ####
## Spring ####
gMod.sig1 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig2 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig3 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig4 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.sig5 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfSp, family = 'nbinom2')
gMod.const <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfSp, family = 'nbinom2')
gMod.no <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfSp, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.nonacy.sp.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)

## Fall ####
gMod.sig1 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +ag_sig1 + impermeable_sig1 + weedy_sig1 + water_sig1 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1) + log(Geocoris+1) + log(Ichneumonoidea+1), # nested random effects not fitted
                     data = dfFa, family = 'nbinom2')
gMod.sig2 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +ag_sig2 + impermeable_sig2 + weedy_sig2 + water_sig2 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig3 <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                       alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +ag_sig3 + impermeable_sig3 + weedy_sig3 + water_sig3 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig4 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +ag_sig4 + impermeable_sig4 + weedy_sig4 + water_sig4 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.sig5 <- glmmTMB(NonAcy ~ Treatment +  wateringMethod + # non-landcover effects
                       alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +ag_sig5 + impermeable_sig5 + weedy_sig5 + water_sig5 + # landcover effects
                       + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                       (1|Site/Field), # nested random effects
                     data = dfFa, family = 'nbinom2')
gMod.const <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                        alfalfa_const + naturalArid_const + dirt_const +ag_const + impermeable_const + weedy_const + water_const + # landcover effects
                        + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                        (1|Site/Field), # nested random effects
                      data = dfFa, family = 'nbinom2')
gMod.no <- glmmTMB(NonAcy ~ Treatment + wateringMethod + # non-landcover effects
                     alfalfa_no + naturalArid_no + dirt_no +ag_no + impermeable_no + weedy_no + water_no + # landcover effects
                     + log(Anthocoridae+1) + log(Arachnida+1) + log(Coccinellidae+1)+ log(Geocoris+1) + log(Ichneumonoidea+1) + # predator effects
                     (1|Site/Field), # nested random effects
                   data = dfFa, family = 'nbinom2')

# dredging
sig1.dredge <- dredge(gMod.sig1, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig2.dredge <- dredge(gMod.sig2, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig3.dredge <- dredge(gMod.sig3, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig4.dredge <- dredge(gMod.sig4, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
sig5.dredge <- dredge(gMod.sig5, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
const.dredge <- dredge(gMod.const, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)
no.dredge <- dredge(gMod.no, m.lim = c(0, 4), fixed = 'cond(Treatment)', trace = 2, cluster = clust)

# rbind and recalc AIC
tab.nb.nonacy.fa.scaled <- rbind(sig1.dredge,sig2.dredge,sig3.dredge,sig4.dredge,sig5.dredge,const.dredge,no.dredge)
