# explores the top down effects of predators. Generates plots.

mDat <- subplotDataRaw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
# Anthocoridae (spring)
anth.eff <- glmmTMB(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                    data = mDat,
                    family = 'nbinom2')
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lmer vers
anth.eff <- lmer(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                 data = mDat)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lm vers
anth.eff <- lm(log(AllAph+1)~ Treatment + log(Anthocoridae+1),
               data = mDat)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
mDat <- subplotDataRaw %>%
  left_join(diffData_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffAnthocoridae > 0)
anth.eff <- glmmTMB(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                    data = mDat,
                    family = 'nbinom2')
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lmer vers
anth.eff <- lmer(AllAph~ Treatment + log(Anthocoridae+1) + (1|Site:Field),
                 data = mDat)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))
# try lm vers
anth.eff <- lm(log(AllAph+1)~ Treatment + log(Anthocoridae+1),
               data = mDat)
summary(anth.eff) # no apparent effect of sham
plot(allEffects(anth.eff))


# Arachnida (spring)
mDat <- subplotDataRaw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
ara.eff <- glmmTMB(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                   data = mDat,
                   family = 'nbinom2')
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lmer vers
ara.eff <- lmer(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                data = mDat)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lm vers
ara.eff <- lm(log(AllAph+1)~ Treatment + log(Arachnida+1),
              data = mDat)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
mDat <- subplotDataRaw %>%
  left_join(diffData_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffArachnida > 0)
ara.eff <- glmmTMB(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                   data = mDat,
                   family = 'nbinom2')
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lmer vers
ara.eff <- lmer(AllAph~ Treatment + log(Arachnida+1) + (1|Site:Field),
                data = mDat)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))
# try lm vers
ara.eff <- lm(log(AllAph+1)~ Treatment + log(Arachnida+1),
              data = mDat)
summary(ara.eff) # no apparent effect of sham
plot(allEffects(ara.eff))


# Coccinellidae (spring)
mDat <- subplotDataRaw %>% filter(Season == 'Spring',
                                  Treatment != 'Pre-')
cocc.eff <- glmmTMB(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                    data = mDat,
                    family = 'nbinom2')
summary(cocc.eff) # no apparent effect of sham
plot(allEffects(cocc.eff))
# try lmer vers
cocc.eff <- lmer(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                 data = mDat)
summary(cocc.eff) # no apparent effect of sham
plot(allEffects(cocc.eff))
# try lm vers
cocc.eff <- lm(log(AllAph+1)~ Treatment + log(Coccinellidae+1),
               data = mDat)
summary(cocc.eff) # no apparent effect of sham
plot(allEffects(cocc.eff))
# literally no effect of sham in any of these models
# constrain data to when sham diff is +0
### THIS is when the coccinellidae effects show up ####
mDat <- subplotDataRaw %>%
  left_join(diffData_wide) %>%
  filter(Season == 'Spring',
         Treatment != 'Pre-',
         diffCoccinellidae > 0)
cocc.eff <- glmmTMB(AllAph~ Treatment + log(Coccinellidae+1) + (1|Site:Field),
                    data = mDat,
                    family = 'nbinom2')
summary(cocc.eff) # -TreatmentSham**
plot(allEffects(cocc.eff))
plot(simulateResiduals(cocc.eff))

## try fall effect of ichneumonoida
mDat.fa <- subplotDataRaw %>%
  left_join(diffData_wide) %>%
  filter(Season == 'Fall',
         Treatment != 'Pre-',
         diffIchneumonoidea > 40) # no treatment effect, regardless of how strict you are here
ich.eff <- glmmTMB(AllAph~ Treatment + log(Ichneumonoidea+1) + (1|Site:Field),
                   data = mDat,
                   family = 'nbinom2')
summary(ich.eff) # no apparent effect of sham
plot(allEffects(ich.eff))
plot(simulateResiduals(ich.eff))




## figure for paper (WIP) ####
anth.tidy <- tidy(anth.eff) %>% filter(term == 'TreatmentSham')
ara.tidy <- tidy(ara.eff) %>% filter(term == 'TreatmentSham')
cocc.tidy <- tidy(cocc.eff) %>% filter(term == 'TreatmentSham') %>%
  select(term:p.value)
ich.tidy <- tidy(ich.eff) %>% filter(term == 'TreatmentSham') %>%
  select(term:p.value)
ncol(anth.tidy)
ncol(ara.tidy)
ncol(cocc.tidy)
ncol(ich.tidy)
tidy.mods <- rbind(anth.tidy, ara.tidy, cocc.tidy, ich.tidy)
tidy.mods$model <- c("Anthocoridae", "Arachnida", "Coccinellidae", "Ichneumonoidea")

dwplot(tidy.mods) %>%
  relabel_predictors(c(TreatmentSham = "Biocontrol effect")) +
  geom_vline(xintercept = 0, color = 'red') +
  theme_classic()
