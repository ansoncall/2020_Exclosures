library(tidyverse)
library(lmerTest)
library(MuMIn)
library(effects)
library(ggeffects)
library(ggfortify)
# import data ####
subplotData <- read_csv('tidy_data/subplotData.csv')

# wrangle ####
# this version currently not used, but this is the original way of
# pooling within fields. This does not work for glms, as AICc cannot be
# calculated when the response contains non-integer data.
fieldData <- subplotData %>%
  filter(Treatment != "Pre-") %>%
  group_by(Site, Field, Season) %>%
  summarize(across(where(is.numeric), .fns = mean)) %>%
  # remove "fixed" data
  select(-contains('_fix')) %>%
  ungroup %>%
  # scale landcover vars
  mutate(across(contains('_'), ~as.numeric(scale(.x)))) %>%
  filter(Season == 'Spring')

# need integers for glm, must use sum(Coccinellidae) instead of
# log(Coccinellidae). this slightly changes the model for gaussian lms as well.
fieldData1 <-subplotData %>%
  filter(Treatment != 'Pre-') %>%
  group_by(Site, Field, Season) %>%
  summarize(Coccinellidae = sum(Coccinellidae),
            wateringMethod = first(wateringMethod)) %>%
  select(Site, Field, Season, Coccinellidae, wateringMethod)

fieldData2 <- subplotData %>%
  filter(Treatment != 'Pre-') %>%
  group_by(Site, Field, Season) %>%
  summarize(across(where(is.numeric), .fns = mean)) %>%
  # remove "fixed" data
  select(-contains('_fix')) %>%
  ungroup %>%
  # scale landcover vars
  mutate(across(contains('_'), ~as.numeric(scale(.x)))) %>%
  select(Site, Field, Season, contains('_'))

fieldData3 <- left_join(fieldData1, fieldData2) %>%
  filter(Season == 'Spring')


# View(fieldData)
# names(fieldData)


# fit global models ####
## make LMER global mods ####
globalLMERsig1 <- lmer(log(Coccinellidae + 1) ~
                         alfalfa_sig1 +
                         naturalArid_sig1 +
                         dirt_sig1 +
                         ag_sig1 +
                         impermeable_sig1 +
                         weedy_sig1 +
                         wet_sig1 +
                         water_sig1 +
                         wateringMethod +
                         (1|Site),
                       data = fieldData3,
                       REML = FALSE,
                       na.action = 'na.fail')
globalLMERsig2 <- lmer(log(Coccinellidae + 1) ~
                         alfalfa_sig2 +
                         naturalArid_sig2 +
                         dirt_sig2 +
                         ag_sig2 +
                         impermeable_sig2 +
                         weedy_sig2 +
                         wet_sig2 +
                         water_sig2 +
                         wateringMethod +
                         (1|Site),
                       data = fieldData3,
                       REML = FALSE,
                       na.action = 'na.fail')
globalLMERsig3 <- lmer(log(Coccinellidae + 1) ~
                         alfalfa_sig3 +
                         naturalArid_sig3 +
                         dirt_sig3 +
                         ag_sig3 +
                         impermeable_sig3 +
                         weedy_sig3 +
                         wet_sig3 +
                         water_sig3 +
                         wateringMethod +
                         (1|Site),
                       data = fieldData3,
                       REML = FALSE,
                       na.action = 'na.fail')
globalLMERsig4 <- lmer(log(Coccinellidae + 1) ~
                         alfalfa_sig4 +
                         naturalArid_sig4 +
                         dirt_sig4 +
                         ag_sig4 +
                         impermeable_sig4 +
                         weedy_sig4 +
                         wet_sig4 +
                         water_sig4 +
                         wateringMethod +
                         (1|Site),
                       data = fieldData3,
                       REML = FALSE,
                       na.action = 'na.fail')
globalLMERsig5 <- lmer(log(Coccinellidae + 1) ~
                         alfalfa_sig5 +
                         naturalArid_sig5 +
                         dirt_sig5 +
                         ag_sig5 +
                         impermeable_sig5 +
                         weedy_sig5 +
                         wet_sig5 +
                         water_sig5 +
                         wateringMethod +
                         (1|Site),
                       data = fieldData3,
                       REML = FALSE,
                       na.action = 'na.fail')
globalLMERconst <- lmer(log(Coccinellidae + 1) ~
                          alfalfa_const +
                          naturalArid_const +
                          dirt_const +
                          ag_const +
                          impermeable_const +
                          weedy_const +
                          wet_const +
                          water_const +
                          wateringMethod +
                          (1|Site),
                        data = fieldData3,
                        REML = FALSE,
                        na.action = 'na.fail')
globalLMERno <- lmer(log(Coccinellidae + 1) ~
                       alfalfa_no +
                       naturalArid_no +
                       dirt_no +
                       ag_no +
                       impermeable_no +
                       weedy_no +
                       wet_no +
                       water_no +
                       wateringMethod +
                       (1|Site),
                     data = fieldData3,
                     REML = FALSE,
                     na.action = 'na.fail')

## make LM global mods ####
globalLMsig1 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig1 +
                     naturalArid_sig1 +
                     dirt_sig1 +
                     ag_sig1 +
                     impermeable_sig1 +
                     weedy_sig1 +
                     wet_sig1 +
                     water_sig1 +
                     wateringMethod,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMsig2 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig2 +
                     naturalArid_sig2 +
                     dirt_sig2 +
                     ag_sig2 +
                     impermeable_sig2 +
                     weedy_sig2 +
                     wet_sig2 +
                     water_sig2 +
                     wateringMethod,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMsig3 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig3 +
                     naturalArid_sig3 +
                     dirt_sig3 +
                     ag_sig3 +
                     impermeable_sig3 +
                     weedy_sig3 +
                     wet_sig3 +
                     water_sig3,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMsig3 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig3 +
                     naturalArid_sig3 +
                     dirt_sig3 +
                     ag_sig3 +
                     impermeable_sig3 +
                     weedy_sig3 +
                     wet_sig3 +
                     water_sig3,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMsig4 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig4 +
                     naturalArid_sig4 +
                     dirt_sig4 +
                     ag_sig4 +
                     impermeable_sig4 +
                     weedy_sig4 +
                     wet_sig4 +
                     water_sig4 +
                     wateringMethod,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMsig5 <- lm(log(Coccinellidae + 1) ~
                     alfalfa_sig5 +
                     naturalArid_sig5 +
                     dirt_sig5 +
                     ag_sig5 +
                     impermeable_sig5 +
                     weedy_sig5 +
                     wet_sig5 +
                     water_sig5 +
                     wateringMethod,
                   data = fieldData3,
                   na.action = 'na.fail')
globalLMconst <- lm(log(Coccinellidae + 1) ~
                      alfalfa_const +
                      naturalArid_const +
                      dirt_const +
                      ag_const +
                      impermeable_const +
                      weedy_const +
                      wet_const +
                      water_const +
                      wateringMethod,
                    data = fieldData3,
                    na.action = 'na.fail')
globalLMno <- lm(log(Coccinellidae + 1) ~
                   alfalfa_no +
                   naturalArid_no +
                   dirt_no +
                   ag_no +
                   impermeable_no +
                   weedy_no +
                   wet_no +
                   water_no +
                   wateringMethod,
                 data = fieldData3,
                 na.action = 'na.fail')

## make GLM global mods ####
globalGLMsig1 <- glm(Coccinellidae ~
                       alfalfa_sig1 +
                       naturalArid_sig1 +
                       dirt_sig1 +
                       ag_sig1 +
                       impermeable_sig1 +
                       weedy_sig1 +
                       wet_sig1 +
                       water_sig1 +
                       wateringMethod,
                     data = fieldData3,
                     family = 'poisson',
                     na.action = 'na.fail')
globalGLMsig2 <- glm(Coccinellidae ~
                       alfalfa_sig2 +
                       naturalArid_sig2 +
                       dirt_sig2 +
                       ag_sig2 +
                       impermeable_sig2 +
                       weedy_sig2 +
                       wet_sig2 +
                       water_sig2 +
                       wateringMethod,
                     data = fieldData3,
                     family = 'poisson',
                     na.action = 'na.fail')
globalGLMsig3 <- glm(Coccinellidae ~
                       alfalfa_sig3 +
                       naturalArid_sig3 +
                       dirt_sig3 +
                       ag_sig3 +
                       impermeable_sig3 +
                       weedy_sig3 +
                       wet_sig3 +
                       water_sig3 +
                       wateringMethod,
                     data = fieldData3,
                     family = 'poisson',
                     na.action = 'na.fail')
globalGLMsig4 <- glm(Coccinellidae ~
                       alfalfa_sig4 +
                       naturalArid_sig4 +
                       dirt_sig4 +
                       ag_sig4 +
                       impermeable_sig4 +
                       weedy_sig4 +
                       wet_sig4 +
                       water_sig4 +
                       wateringMethod,
                     data = fieldData3,
                     family = 'poisson',
                     na.action = 'na.fail')
globalGLMsig5 <- glm(Coccinellidae ~
                       alfalfa_sig5 +
                       naturalArid_sig5 +
                       dirt_sig5 +
                       ag_sig5 +
                       impermeable_sig5 +
                       weedy_sig5 +
                       wet_sig5 +
                       water_sig5 +
                       wateringMethod,
                     data = fieldData3,
                     family = 'poisson',
                     na.action = 'na.fail')
globalGLMconst <- glm(Coccinellidae ~
                        alfalfa_const +
                        naturalArid_const +
                        dirt_const +
                        ag_const +
                        impermeable_const +
                        weedy_const +
                        wet_const +
                        water_const +
                        wateringMethod,
                      data = fieldData3,
                      family = 'poisson',
                      na.action = 'na.fail')
globalGLMno <- glm(Coccinellidae ~
                     alfalfa_no +
                     naturalArid_no +
                     dirt_no +
                     ag_no +
                     impermeable_no +
                     weedy_no +
                     wet_no +
                     water_no +
                     wateringMethod,
                   data = fieldData3,
                   family = 'poisson',
                   na.action = 'na.fail')

# dredge models ####
lmerSig1Dredge <- dredge(globalLMERsig1, m.lim = c(0, 2), trace = 2)
lmerSig2Dredge <- dredge(globalLMERsig2, m.lim = c(0, 2), trace = 2)
lmerSig3Dredge <- dredge(globalLMERsig3, m.lim = c(0, 2), trace = 2)
lmerSig4Dredge <- dredge(globalLMERsig4, m.lim = c(0, 2), trace = 2)
lmerSig5Dredge <- dredge(globalLMERsig5, m.lim = c(0, 2), trace = 2)
lmerConstDredge <- dredge(globalLMERconst, m.lim = c(0, 2), trace = 2)
lmerNoDredge <- dredge(globalLMERno, m.lim = c(0, 2), trace = 2)

lmSig1Dredge <- dredge(globalLMsig1, m.lim = c(0, 2), trace = 2)
lmSig2Dredge <- dredge(globalLMsig2, m.lim = c(0, 2), trace = 2)
lmSig3Dredge <- dredge(globalLMsig3, m.lim = c(0, 2), trace = 2)
lmSig4Dredge <- dredge(globalLMsig4, m.lim = c(0, 2), trace = 2)
lmSig5Dredge <- dredge(globalLMsig5, m.lim = c(0, 2), trace = 2)
lmConstDredge <- dredge(globalLMconst, m.lim = c(0, 2), trace = 2)
lmNoDredge <- dredge(globalLMno, m.lim = c(0, 2), trace = 2)

glmSig1Dredge <- dredge(globalGLMsig1, m.lim = c(0, 2), trace = 2)
glmSig2Dredge <- dredge(globalGLMsig2, m.lim = c(0, 2), trace = 2)
glmSig3Dredge <- dredge(globalGLMsig3, m.lim = c(0, 2), trace = 2)
glmSig4Dredge <- dredge(globalGLMsig4, m.lim = c(0, 2), trace = 2)
glmSig5Dredge <- dredge(globalGLMsig5, m.lim = c(0, 2), trace = 2)
glmConstDredge <- dredge(globalGLMconst, m.lim = c(0, 2), trace = 2)
glmNoDredge <- dredge(globalGLMno, m.lim = c(0, 2), trace = 2)

# rbind model selection tables ####
# lmer
lmerTab <- rbind(
  lmerSig1Dredge,
  lmerSig2Dredge,
  lmerSig3Dredge,
  lmerSig4Dredge,
  lmerSig5Dredge,
  lmerConstDredge,
  lmerNoDredge
)
# lm
lmTab <- rbind(
  lmSig1Dredge,
  lmSig2Dredge,
  lmSig3Dredge,
  lmSig4Dredge,
  lmSig5Dredge,
  lmConstDredge,
  lmNoDredge
)
# glm
glmTab <- rbind(
  glmSig1Dredge,
  glmSig2Dredge,
  glmSig3Dredge,
  glmSig4Dredge,
  glmSig5Dredge,
  glmConstDredge,
  glmNoDredge
)
# all
allMods <- rbind(
  lmerSig1Dredge,
  lmerSig2Dredge,
  lmerSig3Dredge,
  lmerSig4Dredge,
  lmerSig5Dredge,
  lmerConstDredge,
  lmerNoDredge,
  lmSig1Dredge,
  lmSig2Dredge,
  lmSig3Dredge,
  lmSig4Dredge,
  lmSig5Dredge,
  lmConstDredge,
  lmNoDredge,
  glmSig1Dredge,
  glmSig2Dredge,
  glmSig3Dredge,
  glmSig4Dredge,
  glmSig5Dredge,
  glmConstDredge,
  glmNoDredge
)

View(allMods)


# get top models ####
get.models(lmerTab, 1)
get.models(lmTab, 1)
# these are the same
get.models(glmTab, 1)
test <-glmer(Coccinellidae * 3 ~ dirt_sig1 + weedy_sig1 + (1|Site),
             fieldData,
             'poisson')
summary(test)
# these are also the exact same. no random effect variance

# compare lm/glm ####
# build top models
bestLM <- lm(log(Coccinellidae + 1) ~ dirt_sig1 + weedy_sig1, fieldData3)
bestGLM <- glm(Coccinellidae ~ dirt_sig1 + weedy_sig1,
               fieldData3,
               family = 'poisson')

# inspect lm
autoplot(bestLM)
plot(allEffects(bestLM, residuals = TRUE))
lmPredict <- ggpredict(bestLM, terms = c('weedy_sig1', 'dirt_sig1'))
ggplot(lmPredict) +
  geom_line(aes(x, predicted, color = group)) +
  geom_point(data = fieldData3,
             aes(weedy_sig1, Coccinellidae, fill = dirt_sig1),
             shape = 21,
             cex = 5,
             inherit.aes = FALSE) +
  geom_text(data = fieldData3,
            aes(weedy_sig1, Coccinellidae, label=paste0(Site, Field)),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 0) +
  labs(title = 'lm')
summary(bestLM)

# inspect glm
autoplot(bestGLM)
plot(allEffects(bestGLM, residuals = TRUE))
glmPredict <- ggpredict(bestGLM, terms = c('weedy_sig1', 'dirt_sig1'))
ggplot(glmPredict) +
  geom_line(aes(x, predicted, color = group)) +
  geom_point(data = fieldData3,
             aes(weedy_sig1, Coccinellidae, fill = dirt_sig1),
             shape = 21,
             cex = 5,
             inherit.aes = FALSE) +
  geom_text(data = fieldData3,
            aes(weedy_sig1, Coccinellidae, label=paste0(Site, Field)),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 0) +
  labs(title = 'glm')
summary(bestGLM)

model.sel(bestLM, bestGLM)

fieldData3 %>% View

lm