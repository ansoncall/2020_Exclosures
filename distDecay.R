library(tidyverse)
# decay funcs
Distance <- seq(1,1000)
Vaggressive <- c((1.0 / (1.0 + 2.718282 ^ ((3.5 * Distance / 50) - 4.5))) * 1000, rep(0, 100))
Aggressive <- c((1.0 / (1.0 + 2.718282 ^ ((2.0 * Distance / 50) - 6.0))) * 1000, rep(0, 100))
Moderate <- c((1.0 / (1.0 + 2.718282 ^ ((1.25 * Distance / 50) - 8.0))) * 1000, rep(0, 100))
Weak <- c((1.0 / (1.0 + 2.718282 ^ ((1.0 * Distance / 50) - 10.0))) * 1000, rep(0, 100))
Vweak <- c((1.0 / (1.0 + 2.718282 ^ ((0.75 * Distance / 50) - 11.0))) * 1000, rep(0, 100))
Constant <- c((abs(Distance-1000)), rep(0, 100))
None <- c(rep(1000, 1000), rep(0, 100))
# update dist to match length of other cols
Distance <- seq(1, 1100)
dist_decay_df <- data.frame(Distance, Vaggressive, Aggressive, Moderate, Weak, Vweak, Constant, None) %>%
  pivot_longer(-Distance, names_to = "Decay function", values_to = "value") %>%
  mutate(`Decay function` = fct_relevel(`Decay function`, "Vaggressive", "Aggressive", "Moderate", "Weak", "Vweak", "Constant", "None")) %>%
  mutate(`Decay function` = fct_recode(`Decay function`, `Very aggressive` = "Vaggressive", `Very weak` = "Vweak"))

ggplot(dist_decay_df, aes(Distance, value, color = `Decay function`)) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_grey() +
  labs(y = "Weighting value")
