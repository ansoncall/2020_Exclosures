# Exploration of arthropod data
## Aphid histograms ####

# density plot
subplot_data %>%
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

# histogram plot
# all taxa
taxalist <- c("Acyrthosiphon", "NonAcy","Arachnida", "Anthocoridae",
              "Coccinellidae", "Geocoris", "Ichneumonoidea")
subplot_data %>%
  # lengthen (aphids only)
  pivot_longer(all_of(taxalist),
               names_to = "Taxa",
               values_to = "Mean_Density") %>%
  mutate(Taxa = case_when(Taxa == "NonAcy" ~ "Non-Acyrthosiphon aphid",
                          Taxa != "NonAcy" ~ Taxa)) %>%
  mutate(Taxa = fct_relevel(Taxa, "Acyrthosiphon", "Non-Acyrthosiphon aphid")) %>%
  # log-transform
  mutate(Mean_Density = log(Mean_Density + 1)) %>%
  # relevel factors
  mutate(Taxa = fct_relevel(Taxa, "Acyrthosiphon"),
         Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(x = Mean_Density, fill = Season, color = Season)) +
  geom_histogram(alpha = 0.4,
                 position = position_dodge(0.1)) +
  labs(x = "log(Density + 1)",
       y = "Count") +
  theme_grey(base_size = 10) +
  facet_wrap(~ Taxa, ncol = 1, scales = "free") +
  scale_fill_manual(values = seasons_palette) +
  scale_color_manual(values = seasons_palette) +
  xlim(-1, 10.1) +
  theme(#legend.position = c(0.1, 0.1),
    legend.background = element_rect(linetype = 1, color = NA),
    panel.background = element_rect(fill = NA, color = "black"),
    plot.background = element_rect(fill = "white"),
    # panel.grid.major.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color = "black"),
    strip.background.x = element_rect(fill = "NA", color = "NA"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box.margin =  margin(r = 0.2, l = -40, t = 0))

ggsave("season_histogram.pdf", width = 8.5, height = 20, units = "cm",
       dpi = 600)

# final figure ####
subplot_data %>%
  # lengthen (aphids only)
  pivot_longer(all_of(taxalist),
               names_to = "Taxa",
               values_to = "Mean_Density") %>%
  mutate(Taxa = case_when(Taxa == "NonAcy" ~ "Non-Acyrthosiphon aphid",
                          Taxa != "NonAcy" ~ Taxa)) %>%
  mutate(Taxa = fct_relevel(Taxa, "Acyrthosiphon", "Non-Acyrthosiphon aphid")) %>%
  # log-transform
  mutate(Mean_Density4 = Mean_Density * 4) %>% # area in m2
  mutate(Mean_Density = log(Mean_Density4 + 1)) %>%
  # relevel factors
  mutate(Taxa = fct_relevel(Taxa, "Acyrthosiphon"),
         Season = fct_relevel(Season, "Spring")) %>%
  ggplot(aes(x = Mean_Density, y = rev(Taxa), fill = Season, color = Season)) +
  geom_density_ridges(alpha = 0.4,
                      scale = 1,
                      panel_scaling = TRUE) +
  labs(x = expression(paste("log(Density / ", m^2, " + 1)")),
       y = "Count") +
  theme_grey(base_size = 10) +
  facet_wrap(~ Taxa, ncol = 1, scales = "free") +
  scale_fill_manual(values = seasons_palette) +
  scale_color_manual(values = seasons_palette) +
  scale_y_discrete(expand = expansion(add = c(0.05, 0.6))) +
  xlim(-1, 10.1) +
  theme(#legend.position = c(0.1, 0.1),
    legend.background = element_rect(linetype = 1, color = NA),
    panel.background = element_rect(fill = NA, color = "black"),
    plot.background = element_rect(fill = "white"),
    # panel.grid.major.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background.x = element_rect(fill = "NA", color = "NA"),
    legend.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.margin =  margin(r = 0, l = 0, t = 0, b = 0),
    strip.text = element_text(margin = margin(b = 1)))

ggsave("season_density_plot.pdf", width = 8.5, height = 20, units = "cm",
       dpi = 600)
## all taxa - density version


aphlist <- c("Acyrthosiphon", "Aphis", "Therioaphis")
subplot_data %>%
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
subplot_data %>%
  # summarize seasons indepentently
  group_by(Season) %>%
  # create summary cols
  summarize(median_AllAph = median(AllAph),
            median_Acyrthosiphon = median(Acyrthosiphon),
            median_NonAcy = median(NonAcy),
            mean_AllAph = mean(AllAph),
            sd_AllAph = sd(AllAph),
            se_AllAph = sd(AllAph)/sqrt(189),
            se_NonAcy = sd(NonAcy)/sqrt(189),
            se_Geoc = sd(Geocoris)/sqrt(189),
            se_Ara = sd(Arachnida)/sqrt(189),
            mean_Acyrthosiphon = mean(Acyrthosiphon),
            sd_Acyrthosiphon = sd(Acyrthosiphon),
            se_Acyrthosiphon = sd(Acyrthosiphon)/sqrt(189),
            mean_NonAcy = mean(NonAcy),
            sd_NonAcy = sd(NonAcy),
            median_Anth = median(Anthocoridae),
            median_Ara = median(Arachnida),
            median_Cocc = median(Coccinellidae),
            median_Geoc = median(Geocoris),
            median_Ich = median(Ichneumonoidea),
            median_Nab = median(Nabis),
            mean_Anth = mean(Anthocoridae),
            sd_Anth = sd(Anthocoridae),
            mean_Ara = mean(Arachnida),
            sd_Ara = sd(Arachnida),
            mean_Cocc = mean(Coccinellidae),
            sd_Cocc = sd(Coccinellidae),
            mean_Geoc = mean(Geocoris),
            sd_Geoc = sd(Geocoris),
            mean_Ich = mean(Ichneumonoidea),
            sd_Ich = sd(Ichneumonoidea),
            mean_Nab = mean(Nabis)
  ) %>%
  # move all stats into a single col
  pivot_longer(starts_with("mean") | starts_with("median") | starts_with("sd") | starts_with("se", ignore.case = F),
               names_to = "Taxon", values_to = "Value") %>%
  # parse "Taxon" col
  separate(Taxon, c("Stat", "Taxon"), "_", remove = FALSE) %>%
  # arrange rows
  arrange(Taxon, Season, Stat) %>%
  # relocate columns
  relocate(Taxon, Season, Stat) %>%
  # x16 to get value in density/m2
  mutate(Value = Value*4) -> tb
tb %>%
  pivot_wider(names_from = Taxon, values_from = Value) %>%
  tab_df(alternate.rows = TRUE)
tb %>% filter(Stat == "mean") %>%
  pivot_wider(values_from = Value, names_from = Season) %>%
  mutate(change = Fall/Spring) %>%
  filter(Taxon %in% c('Anth', 'Ara', 'Cocc', 'Geoc', 'Ich')) %>% summarize(
    spMean = mean(Spring),
    faMean = mean(Fall),
    changeMean = mean(change),
    spSD = sd(Spring),
    sdFa = sd(Fall),
    seSp = spSD/(sqrt(189)),
    seFa = sdFa/(sqrt(189))
  )
  # print table to viewer


# statistical tests
# aphids
lm(log(AllAph + 1) ~ Season, subplot_data) %>% summary
lm(log(AllAph + 1) ~ Season, subplot_data) %>% confint()
lm(log(Acyrthosiphon + 1) ~ Season, subplot_data) %>% summary
lm(log(NonAcy + 1) ~ Season, subplot_data) %>% summary
# predators
lm(log(Anthocoridae + 1) ~ Season, subplot_data) %>% summary
lm(log(Arachnida + 1) ~ Season, subplot_data) %>% summary
lm(log(Coccinellidae + 1) ~ Season, subplot_data) %>% summary
lm(log(Geocoris + 1) ~ Season, subplot_data) %>% summary
lm(log(Ichneumonoidea + 1) ~ Season, subplot_data) %>% summary
lm(log(Nabis + 1) ~ Season, subplot_data) %>% summary
# print p values with bonferonni correction for predators
summary(lm(log(Anthocoridae + 1) ~ Season, subplot_data))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 5)
summary(lm(log(Arachnida + 1) ~ Season, subplot_data))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 5)
summary(lm(log(Coccinellidae + 1) ~ Season, subplot_data))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 5)
summary(lm(log(Geocoris + 1) ~ Season, subplot_data))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 5)
summary(lm(log(Ichneumonoidea + 1) ~ Season,
           subplot_data))$coefficients[2, 4] %>%
  p.adjust("bonferroni", 5)

# calculate relative skewness of aphid density across seasons
# use datawizard package
# all aphids
subplot_data %>%
  filter(Season == "Spring") %>%
  pull(AllAph) %>%
  datawizard::skewness(.)
subplot_data %>%
  filter(Season == "Fall") %>%
  pull(AllAph) %>%
  datawizard::skewness(.)
# Acyrthosiphon aphids
subplot_data %>%
  filter(Season == "Spring") %>%
  pull(Acyrthosiphon) %>%
  datawizard::skewness(.)
subplot_data %>%
  filter(Season == "Fall") %>%
  pull(Acyrthosiphon) %>%
  datawizard::skewness(.)
# Non-Acyrthosiphon aphids
subplot_data %>%
  filter(Season == "Spring") %>%
  pull(NonAcy) %>%
  datawizard::skewness(.)
subplot_data %>%
  filter(Season == "Fall") %>%
  pull(NonAcy) %>%
  datawizard::skewness(.)

# density boxplot, all aphids combined
subplot_data %>%
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
subplot_data %>%
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

