# explores the top down effects of predators. Generates plots.

lavaan_df <- subplot_data_raw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
# Anthocoridae (spring)
anth.eff <- glmmTMB(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                    data = lavaan_df,
                    family = 'nbinom2')
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lmer vers
anth.eff <- lmer(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                 data = lavaan_df)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lm vers
anth.eff <- lm(log(AllAph+1)~ Treatment + log(Anthocoridae+1),
               data = lavaan_df)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
lavaan_df <- subplot_data_raw %>%
  left_join(diff_data_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffAnthocoridae > 0)
anth.eff <- glmmTMB(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                    data = lavaan_df,
                    family = 'nbinom2')
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lmer vers
anth.eff <- lmer(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                 data = lavaan_df)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lm vers
anth.eff <- lm(log(AllAph+1)~ Treatment + log(Anthocoridae+1),
               data = lavaan_df)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))


# Arachnida (spring)
lavaan_df <- subplot_data_raw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
ara.eff <- glmmTMB(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                   data = lavaan_df,
                   family = 'nbinom2')
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lmer vers
ara.eff <- lmer(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                data = lavaan_df)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lm vers
ara.eff <- lm(log(AllAph+1)~ Treatment + log(Arachnida+1),
              data = lavaan_df)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
lavaan_df <- subplot_data_raw %>%
  left_join(diff_data_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffArachnida > 0)
ara.eff <- glmmTMB(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                   data = lavaan_df,
                   family = 'nbinom2')
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lmer vers
ara.eff <- lmer(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                data = lavaan_df)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lm vers
ara.eff <- lm(log(AllAph+1)~ Treatment + log(Arachnida+1),
              data = lavaan_df)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))


# Coccinellidae (spring)
lavaan_df <- subplot_data_raw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
cocc_eff <- glmmTMB(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                    data = lavaan_df,
                    family = 'nbinom2')
summary(cocc_eff) # no apparent effect of sham
plot(allEffects(cocc_eff))
# try lmer vers
cocc_eff <- lmer(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                 data = lavaan_df)
summary(cocc_eff) # no apparent effect of sham
plot(allEffects(cocc_eff))
# try lm vers
cocc_eff <- lm(log(AllAph+1)~ Treatment + log(Coccinellidae+1),
               data = lavaan_df)
summary(cocc_eff) # no apparent effect of sham
plot(allEffects(cocc_eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
### THIS is when the coccinellidae effects show up ####
lavaan_df <- subplot_data_raw %>%
  left_join(diff_data_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffCoccinellidae > 0)
cocc_eff <- glmmTMB(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                    data = lavaan_df,
                    family = 'nbinom2')
summary(cocc_eff) # -TreatmentSham**
plot(allEffects(cocc_eff))
plot(simulateResiduals(cocc_eff))

## try fall effect of ichneumonoida
lavaan_df.fa <- subplot_data_raw %>%
  left_join(diff_data_wide) %>%
  filter(Season == 'Fall',
         Treatment != 'Pre-',
         diffIchneumonoidea > 40)
# no treatment effect, regardless of how strict you are here^

ich.eff <- glmmTMB(AllAph~ Treatment + log(Ichneumonoidea+1) + (1|Site:Field),
                   data = lavaan_df,
                   family = 'nbinom2')
summary(ich.eff) # no apparent effect of sham
plot(allEffects(ich.eff))
plot(simulateResiduals(ich.eff))




## figure for paper (WIP) ####
anth.tidy <- tidy(anth.eff) %>% filter(term == 'TreatmentSham')
ara.tidy <- tidy(ara.eff) %>% filter(term == 'TreatmentSham')
cocc.tidy <- tidy(cocc_eff) %>% filter(term == 'TreatmentSham') %>%
  select(term:p.value)
ich.tidy <- tidy(ich.eff) %>% filter(term == 'TreatmentSham') %>%
  select(term:p.value)
ncol(anth.tidy)
ncol(ara.tidy)
ncol(cocc.tidy)
ncol(ich.tidy)
tidy.mods <- rbind(anth.tidy, ara.tidy, cocc.tidy, ich.tidy)
tidy.mods$model <- c("Anthocoridae",
                     "Arachnida",
                     "Coccinellidae",
                     "Ichneumonoidea")

dwplot(tidy.mods) %>%
  relabel_predictors(c(TreatmentSham = "Biocontrol effect")) +
  geom_vline(xintercept = 0, color = 'red') +
  theme_classic()
