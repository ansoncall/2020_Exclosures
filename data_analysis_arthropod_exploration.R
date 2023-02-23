# Exploration of arthropod data
## Aphid histograms ####

# density plot
subplotData %>%
  # lengthen (aphids only)
  pivot_longer(c(Acyrthosiphon, NonAcy),
               names_to = "Taxa",
               values_to = "Mean_Density") %>%
  # log-transform
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  # relevel factors
  mutate(Taxa = fct_relevel(Taxa, "Acyrthosiphon"),
         Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(y = Taxa, x = Mean_Density, fill = Taxa)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = "Aphids, Log+1 Transformation",
       x = "log(Density + 1)",
       y = "Taxon") +
  theme(legend.position = "none") +
  facet_wrap(~ Season, nrow = 2) +
  xlim(-1, 10.1) +
  scale_y_discrete(labels = c("Pea aphids +\n Blue alfalfa aphids",
                            "Other taxa")) +
  scale_fill_manual(values = c("#548235", "#3b3838"),
                    guide = "none")

aphlist <- c("Acyrthosiphon", "Aphis", "Therioaphis")
subplotData %>%
  # lengthen (predators only)
  pivot_longer(all_of(aphlist),
               names_to = "Taxa",
               values_to = "Mean_Density") %>%
  # log-transform
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(y = Taxa, x = Mean_Density, fill = Season)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = "Aphids, Log+1 Transformation, all seasons",
       subtitle = "Plot-level density (means of subplots within a plot)",
       x = "log(Density + 1)",
       y = "Taxon")

# calculate median and mean aphid density in each season
subplotData %>%
  # summarize seasons indepentently
  group_by(Season) %>%
  # create summary cols
  summarize(median_AllAph = median(AllAph),
            median_Acyrthosiphon = median(Acyrthosiphon),
            median_NonAcy = median(NonAcy),
            mean_AllAph = mean(AllAph),
            mean_Acyrthosiphon = mean(Acyrthosiphon),
            mean_NonAcy = mean(NonAcy),
            median_Anth = median(Anthocoridae),
            median_Ara = median(Arachnida),
            median_Cocc = median(Coccinellidae),
            median_Geoc = median(Geocoris),
            median_Ich = median(Ichneumonoidea),
            median_Nab = median(Nabis),
            mean_Anth = mean(Anthocoridae),
            mean_Ara = mean(Arachnida),
            mean_Cocc = mean(Coccinellidae),
            mean_Geoc = mean(Geocoris),
            mean_Ich = mean(Ichneumonoidea),
            mean_Nab = mean(Nabis)) %>%
  # move all stats into a single col
  pivot_longer(starts_with("mean") | starts_with("median"),
               names_to = "Taxon", values_to = "Value") %>%
  # parse "Taxon" col
  separate(Taxon, c("Stat", "Taxon"), "_", remove = FALSE) %>%
  # arrange rows
  arrange(Taxon, Season, Stat) %>%
  # relocate columns
  relocate(Taxon, Season, Stat) %>%
  # print table to viewer
  tab_df(alternate.rows = TRUE)

# statistical tests
# aphids
lm(log(AllAph + 1) ~ Season, subplotData) %>% anova
lm(log(Acyrthosiphon + 1) ~ Season, subplotData) %>% summary
lm(log(NonAcy + 1) ~ Season, subplotData) %>% summary
# predators
lm(log(Anthocoridae + 1) ~ Season, subplotData) %>% summary
lm(log(Arachnida + 1) ~ Season, subplotData) %>% summary
lm(log(Coccinellidae + 1) ~ Season, subplotData) %>% summary
lm(log(Geocoris + 1) ~ Season, subplotData) %>% summary
lm(log(Ichneumonoidea + 1) ~ Season, subplotData) %>% summary
lm(log(Nabis + 1) ~ Season, subplotData) %>% summary
# print p values with bonferonni correction for predators
summary(lm(log(Anthocoridae + 1) ~ Season, subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)
summary(lm(log(Arachnida + 1) ~ Season, subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)
summary(lm(log(Coccinellidae + 1) ~ Season, subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)
summary(lm(log(Geocoris + 1) ~ Season, subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)
summary(lm(log(Ichneumonoidea + 1) ~ Season,
           subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)
summary(lm(log(Nabis + 1) ~ Season,
           subplotData))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 6)

# calculate relative skewness of aphid density across seasons
# use datawizard package
# all aphids
subplotData %>%
  filter(Season == "Spring") %>%
  pull(AllAph) %>%
  datawizard::skewness(.)
subplotData %>%
  filter(Season == "Fall") %>%
  pull(AllAph) %>%
  datawizard::skewness(.)
# Acyrthosiphon aphids
subplotData %>%
  filter(Season == "Spring") %>%
  pull(Acyrthosiphon) %>%
  datawizard::skewness(.)
subplotData %>%
  filter(Season == "Fall") %>%
  pull(Acyrthosiphon) %>%
  datawizard::skewness(.)
# Non-Acyrthosiphon aphids
subplotData %>%
  filter(Season == "Spring") %>%
  pull(NonAcy) %>%
  datawizard::skewness(.)
subplotData %>%
  filter(Season == "Fall") %>%
  pull(NonAcy) %>%
  datawizard::skewness(.)

# density boxplot, all aphids combined
subplotData %>%
  # log-transform
  mutate(Mean_Density = log(AllAph + 1)) %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(x = Season, y = Mean_Density, fill = Season)) +
  geom_boxplot() +
  labs(title = "Aphids, Log+1 Transformation",
       subtitle = "Pooled across taxa and split by season",
       x = "Season",
       y = "log(Density + 1)") +
  ylim(c(0, 9)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"),
                    guide = "none")



## Predator histogram ####
# define list of predators
predlist <- c("Arachnida", "Anthocoridae", "Nabis", "Coccinellidae",
              "Geocoris", "Ichneumonoidea")
# build density plot
# density plot + boxplot
subplotData %>%
  # lengthen (predators only)
  pivot_longer(all_of(predlist),
               names_to = "Taxa",
               values_to = "Mean_Density") %>%
  # log-transform
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  # relevel factors
  mutate(Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(y = Taxa, x = Mean_Density, fill = Season)) +
  geom_density_ridges(alpha = 0.6) +
  labs(title = "Predators, Log+1 Transformation, all seasons",
       subtitle = "Plot-level density (means of subplots within a plot)",
       x = "log(Density + 1)",
       y = "Taxon") +
  # theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral")
