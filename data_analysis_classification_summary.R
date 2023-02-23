# for summarizing landcover classification
## Effect of "fixing" alfalfa classification on alfalfa areaScore
fieldData %>%
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
fieldData %>%
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

## weedy/wet binning effect ####
# what is the correlation between weedy and wet cover?
fieldData %>%
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
fieldData %>%
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
fieldData %>%
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
ggplot(data = vegPlots %>% filter(type == "Margin"),
       aes(x = site, y = shan, fill = season)) +
  geom_boxplot() +
  labs(title = "Shannon diversity index",
       x = "Site", y = "Diversity") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
# Plant species richness
ggplot(data = vegPlots %>% filter(type == "Margin"),
       aes(x = site, y = rich, fill = season)) +
  geom_boxplot() +
  labs(title = "Plant species richness",
       x = "Site", y = "Richness") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
# Total plant cover
vegPlots %>%
  filter(type == "Margin") %>%
  mutate(total_cover = select(., 12:132)) %>%
  rowSums(na.rm = TRUE) %>%
  ggplot(aes(x = site, y = total_cover, fill = season)) +
  geom_boxplot() +
  labs(title = "Plant cover %",
       x = "Site", y = "%") +
  scale_fill_manual(values = c("#76db91", "#9e3c21"))
### Notes ####
## Not seeing much variation here. May want to show by field.

## Remember, Yerington will always be missing. Need to fix colors, factor order,
## etc.lib
