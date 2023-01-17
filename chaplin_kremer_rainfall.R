# graph of average annual rainfall at locations of studies cited in 2011
# Chaplin-Kremer review

# rainfall data from:
# developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY

# Chaplin-Kremer 2011 paper:
# Chaplin-Kramer, Rebecca, Megan E. O’Rourke, Eleanor J. Blitzer, and Claire
# Kremen. “A Meta-Analysis of Crop Pest and Natural Enemy Response to Landscape
# Complexity.” Ecology Letters 14, no. 9 (2011): 922–32.
# https://doi.org/10.1111/j.1461-0248.2011.01642.x.

# other notes:

# Locations are approximate. Two sources from the review are missing (both
# dissertations): O'Rourke 2010 and Ramos 2008

library(tidyverse)
library(ggridges)

raw <- read_csv('raw_data/ChaplinKremerRainfall.csv')

ck <- raw %>% filter(Name != 'ROI')

roiPrecip <- raw %>% filter(Name == 'ROI') %>% pull(MeanAnnualPrecip_m)*100

ggplot(ck, aes(MeanAnnualPrecip_m*100)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = roiPrecip, color = 'red') +
  ylab('Number of studies') +
  xlab('Mean annual precipitation, cm') +
  annotate("text", x = 0.6, y = 15,
           label = "Current study", color = 'red', angle = 90) +
  theme_classic()
