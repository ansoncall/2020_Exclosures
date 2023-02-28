# this script identifies which predators are attracted to sham treatments via
# linear mixed-effects models.

## Identify sham-attracted preds ####
### Prepare data ####

# this df contains plot-level DIFFERENCES between sham and control
# this must use UNLOGGED data because subtraction on the log scale doesn't make
# sense
diffData_wide <- subplotDataRaw %>% # long-format counts
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
         .before = 'alfalfa_no') %>%
  # remove original Sham and Control cols to retain only differences
  select(-ends_with('_Sham'), -ends_with('_Control'))
# maybe add scaled+summed pred col here?
# NOTE: this did not end up being useful, but I tried it!

# join diffdata to plotdata
plotWDiff <- left_join(plotData, diffData_wide)
# ready to model? no transforms applied yet
# not all predictor vars are uniformly distributed
dotchart(plotWDiff$Arachnida)
dotchart(plotWDiff$Anthocoridae)
dotchart(plotWDiff$Coccinellidae)
dotchart(plotWDiff$Geocoris)
dotchart(plotWDiff$Ichneumonoidea)
dotchart(plotWDiff$AllAph)
#log+1
dotchart(log(plotWDiff$Arachnida+1))
dotchart(log(plotWDiff$Anthocoridae+1))
dotchart(log(plotWDiff$Coccinellidae+1))
dotchart(log(plotWDiff$Geocoris+1))
dotchart(log(plotWDiff$Ichneumonoidea+1))
dotchart(log(plotWDiff$AllAph+1)) # response

dotchart(plotWDiff$diffArachnida)
dotchart(plotWDiff$diffAnthocoridae)
dotchart(plotWDiff$diffCoccinellidae)
dotchart(plotWDiff$diffGeocoris)
dotchart(plotWDiff$diffIchneumonoidea)
dotchart(plotWDiff$diffAllAph)

### Make models ####
# Mixed effects needed?
# linear and linear mixed-effects often give different results.
# must use mixed effects to account for non-independence
lmer(diffArachnida~1+(1|Site:Field), data = plotWDiff) %>% summary() # no
t.test(plotWDiff$diffArachnida) # same
lmer(diffAnthocoridae~1+(1|Site:Field), data = plotWDiff) %>% summary() # yes
t.test(plotWDiff$diffAnthocoridae) # still yes, but not the same!
lmer(diffCoccinellidae~1+(1|Site:Field), data = plotWDiff) %>% summary() # yes
t.test(plotWDiff$diffCoccinellidae) # same
lmer(diffGeocoris~1+(1|Site:Field), data = plotWDiff) %>% summary() # no
t.test(plotWDiff$diffGeocoris) # still no, but not the same!
lmer(diffIchneumonoidea~1+(1|Site:Field), data = plotWDiff) %>% summary() # no
t.test(plotWDiff$diffIchneumonoidea) # still no, but not the same!

# positive coeff means shams ATTRACT
#### Spring ####
plotDiff.sp <- plotWDiff %>% filter(Season == 'Spring')
dAra.mod <- lmer(diffArachnida~1+(1|Site:Field), data = plotDiff.sp)
dAnt.mod <- lmer(diffAnthocoridae~1+(1|Site:Field), data = plotDiff.sp)
dCoc.mod <- lmer(diffCoccinellidae~1+(1|Site:Field), data = plotDiff.sp)
dGeo.mod <- lmer(diffGeocoris~1+(1|Site:Field), data = plotDiff.sp)
dIch.mod <- lmer(diffIchneumonoidea~1+(1|Site:Field), data = plotDiff.sp)
dAph.mod <- lmer(diffAllAph~1+(1|Site:Field), data = plotDiff.sp)

# collect stats from these mods
diffMods <- list(dAra.mod, dAnt.mod, dCoc.mod, dGeo.mod, dIch.mod, dAph.mod)
diffStats.sp <- tibble(Taxon = c('Ara', 'Ant', 'Coc', 'Geo', 'Ich', 'Aph'),
                       Est = c(0),
                       Ste= c(0),
                       P = c(0))
for (i in 1:6){

  diffStats.sp$Est[[i]] <- summary(diffMods[[i]])$coefficients[[1]]
  diffStats.sp$Ste[[i]] <- summary(diffMods[[i]])$coefficients[[2]]
  diffStats.sp$P[[i]] <- summary(diffMods[[i]])$coefficients[[5]]

}

#### Fall ####
plotDiff.fa <- plotWDiff %>% filter(Season == 'Fall')
dAra.mod <- lmer(diffArachnida~1+(1|Site:Field), data = plotDiff.fa)
dAnt.mod <- lmer(diffAnthocoridae~1+(1|Site:Field), data = plotDiff.fa)
dCoc.mod <- lmer(diffCoccinellidae~1+(1|Site:Field), data = plotDiff.fa)
dGeo.mod <- lmer(diffGeocoris~1+(1|Site:Field), data = plotDiff.fa)
dIch.mod <- lmer(diffIchneumonoidea~1+(1|Site:Field), data = plotDiff.fa)
dAph.mod <- lmer(diffAllAph~1+(1|Site:Field), data = plotDiff.fa)

# collect stats from these mods
diffMods <- list(dAra.mod, dAnt.mod, dCoc.mod, dGeo.mod, dIch.mod, dAph.mod)
diffStats.fa <- tibble(Taxon = c('Ara', 'Ant', 'Coc', 'Geo', 'Ich', 'Aph'),
                       Est = c(0),
                       Ste= c(0),
                       P = c(0))
for (i in 1:6){

  diffStats.fa$Est[[i]] <- summary(diffMods[[i]])$coefficients[[1]]
  diffStats.fa$Ste[[i]] <- summary(diffMods[[i]])$coefficients[[2]]
  diffStats.fa$P[[i]] <- summary(diffMods[[i]])$coefficients[[5]]

}
