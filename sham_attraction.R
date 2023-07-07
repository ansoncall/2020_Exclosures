# this script identifies which predators are attracted to sham treatments via
# linear mixed-effects models.

## Identify sham-attracted preds ####
### Prepare data ####

# this df contains plot-level DIFFERENCES between sham and control
# this must use UNLOGGED data because subtraction on the log scale doesn't make
# sense
diff_data_wide <- subplot_data_raw %>% # long-format counts
  filter(Treatment != 'Pre-') %>% # remove 'Pre-' treatments
  # make separate cols for each treatment (AND taxon)
  pivot_wider(names_from = Treatment, values_from = Arachnida:NonAcy) %>%
  # calculate difference between Sham and Control subplots for each taxon
  mutate(diffArachnida = Arachnida_Sham - Arachnida_Control,
         diffNabis = Nabis_Sham - Nabis_Control,
         diffGeocoris = Geocoris_Sham - Geocoris_Control,
         diffAnthocoridae = Anthocoridae_Sham - Anthocoridae_Control,
         diffCoccinellidae = Coccinellidae_Sham - Coccinellidae_Control,
         diffIchneumonoidea = Ichneumonoidea_Sham - Ichneumonoidea_Control,
         diffAcyrthosiphon = Acyrthosiphon_Sham - Acyrthosiphon_Control,
         diffAphis = Aphis_Sham - Aphis_Control,
         diffTherioaphis = Therioaphis_Sham - Therioaphis_Control,
         diffLygus = Therioaphis_Sham - Therioaphis_Control,
         diffThysanoptera = Thysanoptera_Sham - Thysanoptera_Control,
         diffOther = Other_Sham - Other_Control,
         diffAllAph = AllAph_Sham - AllAph_Control,
         diffNonAcy = NonAcy_Sham - NonAcy_Control,
         .before = 'alfalfa_no')
# maybe add scaled+summed pred col here?
# NOTE: this did not end up being useful, but I tried it!

# join diffdata to plot_data
plot_w_diff <- left_join(plot_data, diff_data_wide)
# ready to model? no transforms applied yet
# not all predictor vars are uniformly distributed
dotchart(plot_w_diff$Arachnida)
dotchart(plot_w_diff$Anthocoridae)
dotchart(plot_w_diff$Coccinellidae)
dotchart(plot_w_diff$Geocoris)
dotchart(plot_w_diff$Ichneumonoidea)
dotchart(plot_w_diff$AllAph)
#log+1
dotchart(log(plot_w_diff$Arachnida+1))
dotchart(log(plot_w_diff$Anthocoridae+1))
dotchart(log(plot_w_diff$Coccinellidae+1))
dotchart(log(plot_w_diff$Geocoris+1))
dotchart(log(plot_w_diff$Ichneumonoidea+1))
dotchart(log(plot_w_diff$AllAph+1)) # response

dotchart(plot_w_diff$diffArachnida)
dotchart(plot_w_diff$diffAnthocoridae)
dotchart(plot_w_diff$diffCoccinellidae)
dotchart(plot_w_diff$diffGeocoris)
dotchart(plot_w_diff$diffIchneumonoidea)
dotchart(plot_w_diff$diffAllAph)

### Make models ####
# try simple cocc. model
lm(Coccinellidae~Treatment, data = subplot_data_raw %>%
     left_join(diff_data_wide) %>%
     filter(Season == "Spring",
            Treatment != "Pre-"))
lm(Coccinellidae~Treatment, data = subplot_data_raw %>%
     left_join(diff_data_wide) %>%
     filter(Season == "Spring",
            Treatment != "Pre-",
            diffCoccinellidae > 0))

1.944*4
1.028*4
1.028/1.944

1.571*4
2.429*4

2.429/1.571

subplot_data_raw %>%
  left_join(diff_data_wide) %>%
  filter(Season == "Spring",
         Treatment != "Pre-") %>%
  summarise(mean(Coccinellidae)*4,
            (sd(Coccinellidae)/sqrt(72))*4)
# Mixed effects needed?
# linear and linear mixed-effects often give different results.
# must use mixed effects to account for non-independence
lmer(diffArachnida~1+(1|Site:Field), data = plot_w_diff) %>% summary() # no
t.test(plot_w_diff$diffArachnida) # same
lmer(diffAnthocoridae~1+(1|Site:Field), data = plot_w_diff) %>% summary() # yes
t.test(plot_w_diff$diffAnthocoridae) # still yes, but not the same!
lmer(diffCoccinellidae~1+(1|Site:Field), data = plot_w_diff) %>% summary() # yes
t.test(plot_w_diff$diffCoccinellidae) # same
lmer(diffGeocoris~1+(1|Site:Field), data = plot_w_diff) %>% summary() # no
t.test(plot_w_diff$diffGeocoris) # still no, but not the same!
lmer(diffIchneumonoidea~1+(1|Site:Field), data = plot_w_diff) %>% summary() # no
t.test(plot_w_diff$diffIchneumonoidea) # still no, but not the same!

# positive coeff means shams ATTRACT
#### Spring ####
plotDiff.sp <- plot_w_diff %>% filter(Season == 'Spring')
dAra.mod <- lmer(diffArachnida~1+(1|Site:Field), data = plotDiff.sp)
dAnt.mod <- lmer(diffAnthocoridae~1+(1|Site:Field), data = plotDiff.sp)
dCoc.mod <- lmer(diffCoccinellidae~1+(1|Site:Field), data = plotDiff.sp)
dGeo.mod <- lmer(diffGeocoris~1+(1|Site:Field), data = plotDiff.sp)
d_ich_mod <- lmer(diffIchneumonoidea~1+(1|Site:Field), data = plotDiff.sp)
dAph.mod <- lmer(diffAllAph~1+(1|Site:Field), data = plotDiff.sp)

# collect stats from these mods
diffMods <- list(dAra.mod, dAnt.mod, dCoc.mod, dGeo.mod, d_ich_mod, dAph.mod)
diff_stats_sp <- tibble(Taxon = c('Ara', 'Ant', 'Coc', 'Geo', 'Ich', 'Aph'),
                       Est = c(0),
                       Ste= c(0),
                       P = c(0))
for (i in 1:6){

  diff_stats_sp$Est[[i]] <- summary(diffMods[[i]])$coefficients[[1]]
  diff_stats_sp$Ste[[i]] <- summary(diffMods[[i]])$coefficients[[2]]
  diff_stats_sp$P[[i]] <- summary(diffMods[[i]])$coefficients[[5]]

}

#### Fall ####
plotDiff.fa <- plot_w_diff %>% filter(Season == 'Fall')
dAra.mod <- lmer(diffArachnida~1+(1|Site:Field), data = plotDiff.fa)
dAnt.mod <- lmer(diffAnthocoridae~1+(1|Site:Field), data = plotDiff.fa)
dCoc.mod <- lmer(diffCoccinellidae~1+(1|Site:Field), data = plotDiff.fa)
dGeo.mod <- lmer(diffGeocoris~1+(1|Site:Field), data = plotDiff.fa)
d_ich_mod <- lmer(diffIchneumonoidea~1+(1|Site:Field), data = plotDiff.fa)
dAph.mod <- lmer(diffAllAph~1+(1|Site:Field), data = plotDiff.fa)

# collect stats from these mods
diffMods <- list(dAra.mod, dAnt.mod, dCoc.mod, dGeo.mod, d_ich_mod, dAph.mod)
diff_stats_fa <- tibble(Taxon = c('Ara', 'Ant', 'Coc', 'Geo', 'Ich', 'Aph'),
                       Est = c(0),
                       Ste= c(0),
                       P = c(0))
for (i in 1:6){

  diff_stats_fa$Est[[i]] <- summary(diffMods[[i]])$coefficients[[1]]
  diff_stats_fa$Ste[[i]] <- summary(diffMods[[i]])$coefficients[[2]]
  diff_stats_fa$P[[i]] <- summary(diffMods[[i]])$coefficients[[5]]

}
summary(dCoc.mod)
summary(d_ich_mod)

# *4 to change untis to m2
diff_stats_sp %>%
  mutate(Est = Est*4,
         Ste = Ste*4,
         Pcorr = p.adjust(P, "bonferroni", 6) ) %>%
  tab_df()

diff_stats_fa %>%
  mutate(Est = Est*4,
         Ste = Ste*4,
         Pcorr = p.adjust(P, "bonferroni", 6) ) %>%
  tab_df()

diff_data_wide %>%
  filter(Season == "Spring") %>%
  summarize(
    mSAra = mean(Arachnida_Sham)*4,
    mCAra = mean(Arachnida_Control)*4,
    pcAra = mSAra/mCAra,
    mSAnth = mean(Anthocoridae_Sham)*4,
    mCAnth = mean(Anthocoridae_Control)*4,
    pcAnth = mSAnth/mCAnth,
    mSCoc = mean(Coccinellidae_Sham)*4,
    mCCoc = mean(Coccinellidae_Control)*4,
    pcCoc = mSCoc/mCCoc
            )

diff_data_wide %>%
  filter(Season == "Fall") %>%
      summarize(
        mSIch = mean(Ichneumonoidea_Sham)*4,
        mCIch = mean(Ichneumonoidea_Control)*4,
        pcIch = mSIch/mCIch,
        medSIch = median(Ichneumonoidea_Sham)*4,
        medCIch = median(Ichneumonoidea_Control)*4,
        pcMed = medSIch/medCIch
  )
