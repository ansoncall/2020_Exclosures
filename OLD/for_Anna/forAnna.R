# for anna
library(tidyverse)
d <- read_csv("tidy_data/vegPlots.csv")

d_present <- d %>%
  mutate(across(12:50, ~ ifelse(.x == 0, 0, 1)))

d_na <- d %>%
  mutate(across(12:50, ~ ifelse(.x == 0, NA, .x)))

encounters <- d_present %>%
  group_by(type, season) %>%
  summarize(across(12:50, ~ sum(.x))) %>%
  pivot_longer(3:41, names_to = "species", values_to = "encounters")

cover <- d_na %>%
  group_by(site, plotnum) %>%
  summarize(across(12:50, ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
  pivot_longer(3:41, names_to = "species", values_to = "cover") %>%
  filter(cover != "NaN") %>%
  group_by(species) %>%
  summarize(across(cover, ~ mean(.x)), .groups = "keep")




# which plants were most commonly encountered?

ggplot(encounters,
       aes(x = reorder(species, -encounters), y = encounters, fill = season)) +
  geom_col() +
  facet_wrap(~type) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("How many times did we encounter this species in a plot?")

ggsave("encounters.jpg", width = 10, height = 7)

ggplot(cover, aes(x = reorder(species, -cover), y = cover)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Given a species was found in a plot, what was the mean % cover?")

ggsave("cover.jpg", width = 10, height = 7)
