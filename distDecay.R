library(tidyverse)
# decay functions
distance <- seq(1, 1000)
proximate <- c((1.0 / (1.0 + 2.718282 ^ ((3.5 * distance / 50) - 4.5))) * 1000,
               rep(0, 100))
near <- c((1.0 / (1.0 + 2.718282 ^ ((2.0 * distance / 50) - 6.0))) * 1000,
          rep(0, 100))
intermediate <- c((1.0 / (1.0 + 2.718282 ^ ((1.25 * distance / 50) - 8.0)))
                  * 1000,
                  rep(0, 100))
distant <- c((1.0 / (1.0 + 2.718282 ^ ((1.0 * distance / 50) - 10.0))) * 1000,
             rep(0, 100))
most_distant <- c((1.0 / (1.0 + 2.718282 ^ ((0.75 * distance / 50) - 11.0)))
                  * 1000,
                  rep(0, 100))
constant <- c((abs(distance - 1000)), rep(0, 100))
none <- c(rep(1000, 1000), rep(0, 100))
# update dist to match length of other cols
distance <- seq(1, 1100)
dist_decay_df <- tibble(distance,
                        Proximate = proximate,
                        Near = near,
                        Intermediate = intermediate,
                        Distant = distant,
                        `Most distant` = most_distant,
                        Constant = constant,
                        None = none) %>%
  pivot_longer(-distance, names_to = "Decay function", values_to = "value")  %>%
  mutate(`Decay function` = fct_relevel(`Decay function`,
                                        "Proximate",
                                        "Near",
                                        "Intermediate",
                                        "Distant",
                                        "Most distant",
                                        "Constant",
                                        "None"))

ggplot(dist_decay_df, aes(distance, value, color = `Decay function`)) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_grey() +
  labs(y = "Weighting value")
