# build and dredge binomial models of coccinellidae presence / absence and
# dredge.

df_fa_bin <- dfFa %>%
  mutate(CBinary = case_when(Coccinellidae == 0 ~ 0,
                             Coccinellidae > 0 ~ 1))

bin_gmod_sig1 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig1 + naturalArid_sig1 + dirt_sig1 +
                           ag_sig1 + impermeable_sig1 + weedy_sig1 +
                           water_sig1, # landcover effects
                         # nested random effects not fitted
                         data = df_fa_bin,
                         family = "binomial",
                         na.action = "na.fail")
bin_gmod_sig2 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig2 + naturalArid_sig2 + dirt_sig2 +
                           ag_sig2 + impermeable_sig2 + weedy_sig2 +
                           water_sig2 + # landcover effects
                           (1 | Site / Field), # nested random effects
                         data = df_fa_bin,
                         family = "binomial",
                         na.action = "na.fail") # convergence warning
bin_gmod_sig3 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig3 + naturalArid_sig3 + dirt_sig3 +
                           ag_sig3 + impermeable_sig3 + weedy_sig3 +
                           water_sig3 + # landcover effects
                           (1 | Site / Field), # nested random effects
                         data = df_fa_bin,
                         family = "binomial",
                         na.action = "na.fail")
bin_gmod_sig4 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig4 + naturalArid_sig4 + dirt_sig4 +
                           ag_sig4 + impermeable_sig4 + weedy_sig4 +
                           water_sig4 + # landcover effects
                           (1 | Site / Field), # nested random effects
                         data = df_fa_bin,
                         family = "binomial",
                         na.action = "na.fail")
bin_gmod_sig5 <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                           wateringMethod + # non-landcover effects
                           alfalfa_sig5 + naturalArid_sig5 + dirt_sig5 +
                           ag_sig5 + impermeable_sig5 + weedy_sig5 +
                           water_sig5 + # landcover effects
                           (1 | Site / Field), # nested random effects
                         data = df_fa_bin,
                         family = "binomial",
                         na.action = "na.fail")
bin_gmod_const <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                            wateringMethod + # non-landcover effects
                            alfalfa_const + naturalArid_const + dirt_const +
                            ag_const + impermeable_const + weedy_const +
                            water_const + # landcover effects
                            (1 | Site / Field), # nested random effects
                          data = df_fa_bin,
                          family = "binomial",
                          na.action = "na.fail")
bin_gmod_no <- glmmTMB(CBinary ~ Treatment + log_AllAph +
                         wateringMethod + # non-landcover effects
                         alfalfa_no + naturalArid_no + dirt_no + ag_no +
                         impermeable_no + weedy_no +
                         water_no + # landcover effects
                         (1 | Site / Field), # nested random effects
                       data = df_fa_bin,
                       family = "binomial",
                       na.action = "na.fail")


dr_sig1 <- dredge(bin_gmod_sig1,
                  fixed = "cond(Treatment)",
                  m.lim = c(0, 2),
                  trace = 2)
dr_sig2 <- dredge(bin_gmod_sig2,
                  fixed = "cond(Treatment)",
                  m.lim = c(0, 2),
                  trace = 2)
dr_sig3 <- dredge(bin_gmod_sig3,
                  fixed = "cond(Treatment)",
                  m.lim = c(0, 2),
                  trace = 2)
dr_sig4 <- dredge(bin_gmod_sig4,
                  fixed = "cond(Treatment)",
                  m.lim = c(0, 2),
                  trace = 2)
dr_sig5 <- dredge(bin_gmod_sig5,
                  fixed = "cond(Treatment)",
                  m.lim = c(0, 2),
                  trace = 2)
dr_const <- dredge(bin_gmod_const,
                   fixed = "cond(Treatment)",
                   m.lim = c(0, 2),
                   trace = 2)
dr_no <- dredge(bin_gmod_no,
                fixed = "cond(Treatment)",
                m.lim = c(0, 2),
                trace = 2)

all_bin_mods <- rbind(dr_sig1,
                    dr_sig2,
                    dr_sig3,
                    dr_sig4,
                    dr_sig5,
                    dr_const,
                    dr_no)
