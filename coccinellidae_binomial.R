# build and dredge binomial models of coccinellidae presence/absence and dredge.

dfFaBin <- dfFa %>%
  mutate(CBinary = case_when(Coccinellidae == 0 ~ 0,
                             Coccinellidae > 0 ~ 1))

bin.gMod.sig1 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +
                           ag_sig1 + impermeable_sig1 + weedy_sig1 +
                           water_sig1, # landcover effects
                         # nested random effects not fitted
                         data = dfFaBin,
                         family = 'binomial',
                         na.action = "na.fail")
bin.gMod.sig2 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +
                           ag_sig2 + impermeable_sig2 + weedy_sig2 +
                           water_sig2 + # landcover effects
                           (1|Site/Field), # nested random effects
                         data = dfFaBin,
                         family = 'binomial',
                         na.action = "na.fail") # convergence warning
bin.gMod.sig3 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +
                           ag_sig3 + impermeable_sig3 + weedy_sig3 +
                           water_sig3 + # landcover effects
                           (1|Site/Field), # nested random effects
                         data = dfFaBin,
                         family = 'binomial',
                         na.action = "na.fail")
bin.gMod.sig4 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +
                           ag_sig4 + impermeable_sig4 + weedy_sig4 +
                           water_sig4 + # landcover effects
                           (1|Site/Field), # nested random effects
                         data = dfFaBin,
                         family = 'binomial',
                         na.action = "na.fail")
bin.gMod.sig5 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +
                           ag_sig5 + impermeable_sig5 + weedy_sig5 +
                           water_sig5 + # landcover effects
                           (1|Site/Field), # nested random effects
                         data = dfFaBin,
                         family = 'binomial',
                         na.action = "na.fail")
bin.gMod.const <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                            wateringMethod + # non-landcover effects
                            alfalfa_const + naturalArid_const + dirt_const +
                            ag_const + impermeable_const + weedy_const +
                            water_const + # landcover effects
                            (1|Site/Field), # nested random effects
                          data = dfFaBin, family = 'binomial')
bin.gMod.no <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                         wateringMethod + # non-landcover effects
                         alfalfa_no + naturalArid_no + dirt_no +ag_no +
                         impermeable_no + weedy_no +
                         water_no + # landcover effects
                         (1|Site/Field), # nested random effects
                       data = dfFaBin,
                       family = 'binomial',
                       na.action = "na.fail")


dr.sig1 <- dredge(bin.gMod.sig1,
                  fixed = 'cond(Treatment)',
                  m.lim = c(0,2),
                  trace =2)
dr.sig2 <- dredge(bin.gMod.sig2,
                  fixed = 'cond(Treatment)',
                  m.lim = c(0,2),
                  trace =2)
dr.sig3 <- dredge(bin.gMod.sig3,
                  fixed = 'cond(Treatment)',
                  m.lim = c(0,2),
                  trace =2)
dr.sig4 <- dredge(bin.gMod.sig4,
                  fixed = 'cond(Treatment)',
                  m.lim = c(0,2),
                  trace =2)
dr.sig5 <- dredge(bin.gMod.sig5,
                  fixed = 'cond(Treatment)',
                  m.lim = c(0,2),
                  trace =2)
dr.const <- dredge(bin.gMod.const,
                   fixed = 'cond(Treatment)',
                   m.lim = c(0,2),
                   trace =2)
dr.no <- dredge(bin.gMod.no,
                fixed = 'cond(Treatment)',
                m.lim = c(0,2), trace =2)

allBinMods <- rbind(dr.sig1,
                    dr.sig2,
                    dr.sig3,
                    dr.sig4,
                    dr.sig5,
                    dr.const,
                    dr.no)