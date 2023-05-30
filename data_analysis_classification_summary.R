# for summarizing landcover classification
## Effect of "fixing" alfalfa classification on alfalfa areaScore
field_data %>%
  # only need one season
  filter(Season == "Spring") %>%
  # calculate change in alfalfa area
  mutate(Change_sig1 = alfalfa_fix_sig1 - alfalfa_sig1,
         Change_sig2 = alfalfa_fix_sig2 - alfalfa_sig2,
         Change_sig3 = alfalfa_fix_sig3 - alfalfa_sig3,
         Change_sig4 = alfalfa_fix_sig4 - alfalfa_sig4,
         Change_sig5 = alfalfa_fix_sig5 - alfalfa_sig5,
         Change_const = alfalfa_fix_const - alfalfa_const,
         Change_no = alfalfa_fix_no - alfalfa_no) %>%
  # scale across new cols (must ungroup first)
  ungroup() %>%
  mutate(across(contains("Change"), ~as.numeric(scale(.x)))) %>%
  # pivot longer
  pivot_longer(contains("Change"),
               names_to = "distWeight",
               values_to = "Change") %>%
  # make field id, also trim up values in distWeight column
  mutate(FieldID = paste0(Site, "0", Field),
         distWeight = str_extract(distWeight, "_[:alnum:]+")) %>%
  # reorder_within to make ordered bars within facets of ggplot
  mutate(FieldID = reorder_within(FieldID,
                                  desc(Change),
                                  distWeight,
                                  sep = "")) %>%
  ggplot(aes(FieldID, Change)) +
  geom_col() +
  facet_wrap(~ distWeight, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Effect of "fixing" alfalfa classification on weedyWet areaScore
field_data %>%
  # only need one season
  filter(Season == "Spring", Treatment == "Control") %>%
  # calculate change in alfalfa area
  mutate(Change_sig1 = weedyWet_fix_sig1 - weedyWet_sig1,
         Change_sig2 = weedyWet_fix_sig2 - weedyWet_sig2,
         Change_sig3 = weedyWet_fix_sig3 - weedyWet_sig3,
         Change_sig4 = weedyWet_fix_sig4 - weedyWet_sig4,
         Change_sig5 = weedyWet_fix_sig5 - weedyWet_sig5,
         Change_const = weedyWet_fix_const - weedyWet_const,
         Change_no = weedyWet_fix_no - weedyWet_no) %>%
  # scale across new cols (must ungroup first)
  ungroup() %>%
  mutate(across(contains("Change"), ~as.numeric(scale(.x)))) %>%
  # pivot longer
  pivot_longer(contains("Change"),
               names_to = "distWeight",
               values_to = "Change") %>%
  # make field id, also trim up values in distWeight column
  mutate(FieldID = paste0(Site, "0", Field),
         distWeight = str_extract(distWeight, "_[:alnum:]+")) %>%
  # reorder_within to make ordered bars within facets of ggplot
  mutate(FieldID = reorder_within(FieldID,
                                  desc(Change),
                                  distWeight,
                                  sep = "")) %>%
  ggplot(aes(FieldID, Change)) +
  geom_col() +
  facet_wrap(~ distWeight, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# in summary, manual classification seems to change Yerington the most

## figure: "fixing" effect on yerington 2
field_data %>%
  # only need one season
  filter(Season == "Spring", Treatment == "Control",
         Site == "Yerington", Field == 2) %>%
  # calculate change in alfalfa area
  mutate(Change_sig1 = alfalfa_fix_sig1 - alfalfa_sig1,
         Change_sig2 = alfalfa_fix_sig2 - alfalfa_sig2,
         Change_sig3 = alfalfa_fix_sig3 - alfalfa_sig3,
         Change_sig4 = alfalfa_fix_sig4 - alfalfa_sig4,
         Change_sig5 = alfalfa_fix_sig5 - alfalfa_sig5,
         Change_const = alfalfa_fix_const - alfalfa_const,
         Change_no = alfalfa_fix_no - alfalfa_no) %>%
  rowwise() %>%
  # change areascores to proportions
  mutate(
         across(ends_with("_sig1"), ~ .x/rowSums(across(ends_with("_sig1")))),
         across(ends_with("_sig2"), ~ .x/rowSums(across(ends_with("_sig2")))),
         across(ends_with("_sig3"), ~ .x/rowSums(across(ends_with("_sig3")))),
         across(ends_with("_sig4"), ~ .x/rowSums(across(ends_with("_sig4")))),
         across(ends_with("_sig5"), ~ .x/rowSums(across(ends_with("_sig5")))),
         across(ends_with("_const"), ~ .x/rowSums(across(ends_with("_const")))),
         across(ends_with("_no"), ~ .x/rowSums(across(ends_with("_no")))),
         ) %>%
  ungroup %>%
  select(Site,
         Field,
         contains("alfalfa"),
         contains("weedy"),
         -contains("Perim"),
         -contains("Wet")) %>%
  pivot_longer(contains(c("alfalfa","weedy")),
               names_to = "dist_weight",
               values_to = "proportion") %>%
  separate(dist_weight, c("cover", "dist"), sep = "_(?!f)") %>%
  separate(cover, c("cover", "fix"), sep = "_", fill = "right") %>%
  mutate(fix = replace_na(fix, "regular"),
         dist = fct_relevel(dist,
                            "sig1",
                            "sig2",
                            "sig3",
                            "sig4",
                            "sig5",
                            "const",
                            "no"
                            ),
         dist = fct_recode( dist,
                            `Very aggressive` = "sig1",
                            `Aggressive` = "sig2",
                            `Moderate` = "sig3",
                            `Weak` = "sig4",
                            `Very weak` = "sig5",
                            `Constant` = "const",
                            `None` = "no"),
         fix = fct_relevel(fix, "regular", "fix"),
         fix = fct_recode(fix,
                          `Automatically classified` = "regular",
                          `Manually classified` = "fix"),
         cover = fct_recode(cover,
                            `Alfalfa cover` = "alfalfa",
                            `Weedy cover` = "weedy")) %>%
  filter(cover == "Weedy cover") %>%
    ggplot(aes(dist, proportion*100, fill = fix)) +
    geom_col(
      position = position_dodge(0.9)
      ) +
    # facet_wrap(~ cover, ncol = 1) +
    labs(x = "Distance weighting",
         y = "% weighted area") +
    theme(legend.background = element_rect(linetype = 1, color = NA),
          legend.position = c(0.68, 0.82),
          panel.background = element_rect(fill = NA, color = "black"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 35, hjust = 0.9),
          strip.background.x = element_rect(fill = "NA", color = "NA"),
          legend.title = element_blank())

ggsave("corrected_classifications.pdf", height = 8.5, width = 8.5, units = "cm",
       dpi = 600)

scales::hue_pal()(2)

## ladybug/weedy and fixing
subplot_data_raw %>%
  # spring only
  filter(Season == "Spring") %>%
  # log-transform AllAph col
  rowwise() %>%
  # change areascores to proportions
  mutate(
    across(contains("fix") & ends_with("_sig1"), ~ .x/rowSums(across(contains("fix") & ends_with("_sig1")))),
    across(contains("fix") & ends_with("_sig2"), ~ .x/rowSums(across(contains("fix") & ends_with("_sig2")))),
    across(contains("fix") & ends_with("_sig3"), ~ .x/rowSums(across(contains("fix") & ends_with("_sig3")))),
    across(contains("fix") & ends_with("_sig4"), ~ .x/rowSums(across(contains("fix") & ends_with("_sig4")))),
    across(contains("fix") & ends_with("_sig5"), ~ .x/rowSums(across(contains("fix") & ends_with("_sig5")))),
    across(contains("fix") & ends_with("_const"), ~ .x/rowSums(across(contains("fix") & ends_with("_const")))),
    across(contains("fix") & ends_with("_no"), ~ .x/rowSums(across(contains("fix") & ends_with("_no")))),
  ) %>%
  mutate(
    across(!contains("fix") & ends_with("_sig1"), ~ .x/rowSums(across(!contains("fix") & ends_with("_sig1")))),
    across(!contains("fix") & ends_with("_sig2"), ~ .x/rowSums(across(!contains("fix") & ends_with("_sig2")))),
    across(!contains("fix") & ends_with("_sig3"), ~ .x/rowSums(across(!contains("fix") & ends_with("_sig3")))),
    across(!contains("fix") & ends_with("_sig4"), ~ .x/rowSums(across(!contains("fix") & ends_with("_sig4")))),
    across(!contains("fix") & ends_with("_sig5"), ~ .x/rowSums(across(!contains("fix") & ends_with("_sig5")))),
    across(!contains("fix") & ends_with("_const"), ~ .x/rowSums(across(!contains("fix") & ends_with("_const")))),
    across(!contains("fix") & ends_with("_no"), ~ .x/rowSums(across(!contains("fix") & ends_with("_no")))),
  ) %>%
  ungroup() %>%
  mutate(log_Coccinellidae = log(Coccinellidae + 1)) %>%
  select(log_Coccinellidae, Coccinellidae, weedy_sig1, weedy_fix_sig1) %>%
  pivot_longer(contains("weedy"), names_to = "fix", values_to = "area") %>%
  ggplot(aes(area, log_Coccinellidae, color = fix)) +
  geom_point() +
  geom_smooth(method = "lm")

## weedy/wet binning effect ####
# what is the correlation between weedy and wet cover?
field_data %>%
  # only need one season and trt
  filter(Season == "Spring", Treatment == "Control") %>%
  # focus on weedy vs. wet
  select(contains(c("weedy", "wet"))) %>%
  # drop data from "fixed" classifications, also drop "weedyWet" cols
  select(-contains(c("fix", "weedyWet"))) %>%
  # pivot longer to put all distweights in one col
  pivot_longer(contains(c("weedy", "wet"))) %>%
  # separate klass and distweight in name col
  separate(name, c("klass", "distWeight"), "_") %>%
  # widen
  pivot_wider(names_from = c(klass), values_from = value) %>%
  ggplot(aes(x = weedy, y = wet)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ distWeight, scales = "free") +
  labs(title = 'Correlation between "weedy" and "wet" classes',
       subtitle = "across all distance weights") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

# generally poor correlation between these two classes
# bin these because they are biologically similar, not because they are
# correlated

# wet and water correlation
field_data %>%
  # only need one season and trt
  filter(Season == "Spring", Treatment == "Control") %>%
  # focus on weedy vs. wet
  select(contains(c("water", "wet"))) %>%
  # drop data from "fixed" classifications, also drop "weedyWet" cols
  select(-contains(c("fix"))) %>%
  # pivot longer to put all distweights in one col
  pivot_longer(contains(c("water", "wet"))) %>%
  # separate klass and distweight in name col
  separate(name, c("klass", "distWeight"), "_") %>%
  # widen
  pivot_wider(names_from = c(klass), values_from = value) %>%
  ggplot(aes(x = wet, y = water)) +
  geom_point() +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(~ distWeight, scales = "free") +
  labs(title = 'Correlation between "wet" and "water" classes',
       subtitle = "across all distance weights") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
# moderate correlation here. this breaks down when using weedyWet (7-class)
field_data %>%
  # only need one season and trt
  filter(Season == "Spring", Treatment == "Control") %>%
  # focus on weedy vs. wet
  select(contains(c("water", "weedyWet"))) %>%
  # drop data from "fixed" classifications, also drop "weedyWet" cols
  select(-contains(c("fix"))) %>%
  # pivot longer to put all distweights in one col
  pivot_longer(contains(c("water", "weedyWet"))) %>%
  # separate klass and distweight in name col
  separate(name, c("klass", "distWeight"), "_") %>%
  # widen
  pivot_wider(names_from = c(klass), values_from = value) %>%
  ggplot(aes(x = weedyWet, y = water)) +
  geom_point() +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(~ distWeight, scales = "free") +
  labs(title = 'Correlation between "weedyWet" and "water" classes',
       subtitle = "across all distance weights") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())


## pca ####

# Not wanting to fix this right now. Used to have pca plot for every combination
# of distWeight, nClasses (7 or 8), and fixed vs. non-fixed


# Margin data summary ####
## Across sites ####
# Shannon diversity
ggplot(data = veg_plots %>% filter(type == "Margin"),
       aes(x = site, y = shan, fill = season)) +
  geom_boxplot() +
  labs(title = "Shannon diversity index",
       x = "Site", y = "Diversity") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
# Plant species richness
ggplot(data = veg_plots %>% filter(type == "Margin"),
       aes(x = site, y = rich, fill = season)) +
  geom_boxplot() +
  labs(title = "Plant species richness",
       x = "Site", y = "Richness") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
# Total plant cover
veg_plots %>%
  filter(type == "Margin") %>%
  mutate(total_cover = select(., 12:132) %>%
           rowSums(na.rm = TRUE)) %>%
  ggplot(aes(x = site, y = total_cover, fill = season)) +
  geom_boxplot() +
  labs(title = "Plant cover %",
       x = "Site", y = "%") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
### Notes ####
## Not seeing much variation here. May want to show by field.

## Remember, Yerington will always be missing. Need to fix colors, factor order,
## etc.lib
