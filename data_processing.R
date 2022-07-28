# Title: Data processing for 2020 Exclosures project
# Author: Anson Call
# Date: 2022/02/28
# This script cleans the raw data in preparation for analysis. It only needs to
# be run once.

# load packages ####
library(tidyverse) # R packages for data science
library(magrittr) # for assignment and exposition pipes
library(vegan) # for diversity metrics

# define functions ####
# date-to-season function for vegsurvey data
getSeason <- function(DATES) {
  # define cutoff between spring and fall surveys
  cutoff <- as.Date("2021-6-30",  format = "%Y-%m-%d") # spring end
  # convert dates to proper format
  d <- as.Date(strftime(DATES, format="%Y-%m-%d"))
  # return season
  if_else (d <= cutoff, 'Spring', 'Fall')
}

# import data ####
## spring data ####
spring <-
  read_csv('raw_data/Spring_counts.csv',
           col_types = 'cfffffffddddddddddddc') %>%
  select(-Sorter, -Counter, -Counter_2) %>%
  mutate(AllAph = Acyrthosiphon + Aphis + Therioaphis,
         NonAcy = Aphis + Therioaphis) %>%
  mutate(Season = 'Spring')

## fall data ####
fall <- read_csv('raw_data/Fall_counts.csv',
                 col_types = 'dcffffffdddddddddddd') %>%
  select(-Sorter, -Counter) %>%
  mutate(AllAph = Acyrthosiphon + Aphis + Therioaphis,
         NonAcy = Aphis + Therioaphis) %>%
  mutate(Season = 'Fall')

## veg data ####
vegdata <- read_csv('raw_data/VegSurvey.csv',
                    col_types = 'Dcfffddccc') %>%
  mutate(waypointID = paste(site, waypoint))

## updated plant ids for vegdata ####
new_ids <- read_csv('raw_data/updated_plant_ids.csv') %>%
  select(orig_code, new_code)

## field location data ####
## incl. closest field info for survey points in vegdata
rawlocs <- read_csv('raw_data/vegSurvey_fieldJoin.csv',
                    col_types = 'fffTdfddfcc')

## veg survey plot type data
veg_survey_classes <- read_csv('raw_data/vegSurveyClasses.csv',
                               col_types = 'ffff')

## landcover class data
veg_survey_cover <- read_csv('raw_data/vegData_coverJoin.csv',
                             col_types = 'fffffffffffffffffffffffffcc')

## landcover data ####
# now only 1 csv (supervised classification, 8 classes)
# import csv files
landcover <- read_csv('raw_data/superDove_areaScores.csv')

# check data ####
## check spring and fall data ####
spring %>%
  group_by(Site, Field) %>%
  count()
spring %>%
  group_by(Site, Treatment) %>%
  count()

# Identify error 1
spring %>%
  group_by(Site, Field, Treatment) %>%
  filter(Site=='Lovelock')
#L3P2Full recorded as field 2 instead of field 3

# Resolve error 1: L3P2Full field 2 -> field 3
spring %<>%
  mutate(Field = case_when(Vial=='L3P2Full' ~ '3',
                           Vial!='L3P2Full' ~ as.character(Field)))

# Identify error 2
spring %>%
  group_by(Site, Field, Treatment) %>% filter(Site=='Lovelock') %>% count()
# Lovelock field 3 now overweight
spring %>%
  group_by(Site, Field, Treatment) %>% filter(Site=='Lovelock')
fall %>%
  group_by(Site, Field, Treatment) %>% filter(Site=='Lovelock')
# L3P2Full (and L3P3SHAM) match text pattern of fall data
fall %>%
  group_by(Site, Field, Treatment) %>% filter(Site=='Lovelock') %>% count()
# Lovelock field 3 full underweight
# suggests spring L3P2Full should actually be fall F3P2Full. check data
spring %>%
  filter(Site=='Lovelock', Field==3) %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
spring %>%
  filter(Vial=='L3P2Full') %>%
  select(AllAph) # 465
# L3P2Full aphid count is way higher than other Lovelock field 3 full counts in
# spring data further suggests L3P2Full belongs in the fall data
fall %>%
  filter(Site=='Lovelock', Field==3) %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  geom_abline(slope = 0, intercept = 465)
# L3P2Full aphid count fits in with fall data well

# Resolve error 2: spring L3P2Full -> fall
spring %<>%
  mutate(Season = case_when(Vial=='L3P2Full' ~ 'Fall',
                            Vial!='L3P2Full' ~ as.character(Season)))

# Identify error 3
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
# Lovelock field 3 overweight
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
# spring Lovelock sham overweight by 2, Lovelock control underweight by 1
fall %>%
  group_by(Site, Treatment) %>%
  filter(Site=='Lovelock') %>%
  count
# fall Lovelock sham is underweight by 1
# suggests at least 1 spring Lovelock sham belongs in fall data
spring %>%
  group_by(Site, Field, Treatment) %>%
  filter(Site=='Lovelock', Treatment=='Sham')
fall %>%
  group_by(Site, Field, Treatment) %>%
  filter(Site=='Lovelock', Treatment=='Sham')
# spring L3P3SHAM matches text pattern of fall data
# further suggests L3P3SHAM belongs in the fall data
spring %>%
  filter(Site=='Lovelock', Field==3) %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
spring %>%
  filter(Vial=='L3P3SHAM') %>%
  select(AllAph) # 709
# L3P3SHAM aphid count is way higher than other Lovelock field 3 full counts in
# spring data
fall %>%
  filter(Site=='Lovelock', Field==3) %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  geom_abline(slope = 0, intercept = 709)
# L3P3SHAM aphid count is still high for fall data, but fall data is still
# higher, on average, than spring data. further suggests L3P3SHAM is a better
# fit in the fall data

# Resolve error 3: spring L3P3SHAM -> fall
spring %<>%
  mutate(Season = case_when(Vial=='L3P3SHAM' ~ 'Fall',
                            Vial!='L3P3SHAM' ~ as.character(Season)))

# Identify error 4
spring %>%
  group_by(Site, Field) %>%
  filter(Season=='Spring') %>%
  count
# Yerington field 2?1 measurements are probably field 1, since Yerington field 1
# is underweight by 2
spring %>%
  group_by(Site, Field, Treatment) %>%
  filter(Season=='Spring', Site=='Yerington') %>%
  count
# Yerington field 2?1 treatments match the treatments that are underweight in
# Yerington field 1 further suggests that 2?1 data could belong in field 1
spring %>%
  filter(Site=='Yerington', Field=='1' | Field=='2' | Field=='2?1') %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Field)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# YF2(YF1?)P2SHAM is actually two rows (!)
spring %>%
  filter(Vial=='YF2(YF1?)P2SHAM')
# Which one is the true field 2 sham?
spring %>%
  filter(Site=='Yerington', Field=='1' | Field=='2' | Field=='2?1') %>%
  mutate(VialID = paste(Vial, Field)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# field 2 generally has higher aphid counts, so we should expect YF2(YF1?)P2SHAM
# 2?1 to be the one from field 1, since it has a lower count

# Correct error 4: YF2(YF1?)P2SHAM field 2?1 -> field 1
spring %<>%
  mutate(Field = case_when(Vial=='YF2(YF1?)P2SHAM' & Field=='2?1' ~ '1',
                           Vial!='YF2(YF1?)P2SHAM' | Field!='2?1' ~
                             as.character(Field)))

# Identify error 5
spring %>%
  group_by(Site, Field) %>%
  filter(Season=='Spring') %>%
  count
# 1 Yerington field 2?1 measurement remaining, Yerington field 1 still
# underweight by 1
spring %>%
  group_by(Site, Field, Treatment) %>%
  filter(Season=='Spring',
         Site=='Yerington',
         Field=='1' | Field=='2' | Field=='2?1') %>%
  count
# remaining 2?1 measurement is a full measurement, matching what is missing from
# Yerington field 1
spring %>%
  filter(Site=='Yerington', Field=='1' | Field=='2' | Field=='2?1') %>%
  mutate(VialID = paste(Vial, Field)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# remaining 2?1 measurement has low aphid density, matching other measurements
# from field 1 further suggests that the remaining 2?1 measurement is from field
# 1

# Correct error 5: YF2(YF1?)P2FULL field 2?1 -> field 1
spring %<>%
  mutate(Field = case_when(Vial=='YF2(YF1?)P2FULL' & Field=='2?1' ~ '1',
                           Vial!='YF2(YF1?)P2FULL' | Field!='2?1' ~
                             as.character(Field)))

# Identify error 6
spring %>%
  group_by(Site, Field) %>%
  filter(Season=='Spring') %>%
  count
# missing 2 measurements from Fallon field 4
spring %>%
  group_by(Field, Plot) %>%
  filter(Season=='Spring', Site=='Fallon') %>%
  count
# both missing measurements are from Fallon field 4, plot 3
spring %>%
  filter(Season=='Spring', Site=='Fallon', Field=='4', Plot=='3') %>%
  arrange(Treatment)
# Fallon field 4, plot 3 full and sham are missing
# Are they in the fall dataset?
fall %>%
  group_by(Site, Field) %>%
  count
# fall Fallon 4 contains 2 extra measurements
fall %>%
  group_by(Site, Field) %>%
  filter(Site=='Fallon', Field=='4') %>%
  arrange(Treatment)
# hard to tell which measurements are out of place because text patterns in Vial
# names are already standardized
# duplicates of field 4 plot 3 full (one of these may be plot 2 or fall)
# duplicates of field 4 plot 3 sham (one of these may be fall)
# duplicates of plot 2 pre- (one of these may be full?)
# examine fall data to see if there are any outliers that suggest a measurement
# belongs in the spring data
fall %>%
  filter(Site=='Fallon', Field=='4') %>%
  mutate(VialID = paste(Vial, AllAph)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('Pre-' = 'red', 'Full' = 'green',
                               'Control' = 'blue', 'Sham' = 'purple'))
# 1 FALF4P3 sham measurement seems low
# does it fit better with the spring data?
spring %>%
  filter(Season=='Spring', Site=='Fallon', Field=='4') %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  geom_abline(slope = 0, intercept = 18) +
  scale_fill_manual(values = c('Pre-' = 'red', 'Full' = 'green',
                               'Control' = 'blue', 'Sham' = 'purple'))
# This measurement seems low compared to other field 4 sham measurements in both
# spring and fall
# What about plot 3 measurements across spring and fall?
fall %>%
  filter(Site=='Fallon', Field=='4', Plot=='3') %>%
  mutate(VialID = paste(Vial, AllAph)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('Pre-' = 'red', 'Full' = 'green',
                               'Control' = 'blue', 'Sham' = 'purple'))
spring %>% filter(Season=='Spring', Site=='Fallon', Field=='4', Plot=='3') %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity') +
  geom_abline(slope = 0, intercept = 77) +
  geom_abline(slope = 0, intercept = 174) +
  scale_fill_manual(values = c('Pre-' = 'red', 'Full' = 'green',
                               'Control' = 'blue', 'Sham' = 'purple'))
# spring has the higher densities for pre- and sham treatments, so it seems
# slightly more likely that the very low FALF4P3 sham measurement and lower of
# the two FALF4P3 full measurements should be retained in the fall, and the two
# higher FALF4P3 measurements should be moved to the spring data

# admittedly, this is a best guess, and could be incorrect - the safest bet
# would be to throw out all measurements that can't be reliably assigned to
# spring or fall data

# however, the overall difference between spring and fall seems negligible, so
# we may end up pooling this data anyway i'll leave it in for now, making
# corrections based on maximum likelihood

# also, there is a question of which measurement belongs in the spring Fallon
# field 4 plot 2 group:

# Correct error 6: fall Fallon field 4, plot 3 sham and full measurements with
# highest aphid density -> spring
fall %<>%
  mutate(Season = case_when(Vial=='FALF4P3 full' & AllAph == 106 ~ 'Spring',
                            Vial!='FALF4P3 full' | AllAph != 106 ~
                              as.character(Season))) %>%
  mutate(Season = case_when(Vial=='FALF4P3 sham' & AllAph == 92 ~ 'Spring',
                            Vial!='FALF4P3 sham' | AllAph != 92 ~
                              as.character(Season)))

# Identify error 7
fall %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Fall') %>%
  count
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Fall') %>%
  count
## fall is missing a Fallon full measurement, and contains one extra Fallon pre-
## measurement
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
fall %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
# spring is missing a Lovelock control measurement, and contains one extra
# Lovelock sham measurement
# Resolve fall data first
fall %>%
  group_by(Site, Field, Treatment) %>%
  filter(Site=='Fallon', Season=='Fall') %>%
  count
# perhaps one Fallon field 4 pre- measurement should move to Fallon field 4
# full?
fall %>%
  group_by(Site, Field, Plot, Treatment) %>%
  filter(Site=='Fallon', Field=='4', Season=='Fall') %>%
  count
# looks like one Fallon field 4 plot 2 pre- measurement should move to Fallon
# field 4 plot 2 full
fall %>%
  filter(Site=='Fallon', Field=='4', Season=='Fall') %>%
  mutate(VialID = paste(Vial, AllAph)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
# I'm going to assume that the lower of the two pre- measurements should be
# changed to full, since pre- measurements generally have the highest aphid
# counts because they represent collections from a larger area. Higher aphid
# counts from pre- treatments aren't very obvious at Fallon Field 4, but
# experiment-wide this is the general trend.

# Resolve error 7: Move the lower Fallon field 4, plot 2 pre- measurement to
# Fallon field 4, plot 2 full
fall %<>%
  mutate(Treatment = case_when(Vial=='FALF4P2' & AllAph == 77 ~ 'Full',
                               Vial!='FALF4P2' | AllAph != 77 ~
                                 as.character(Treatment)))
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Fall') %>%
  count

# Identify error 8
spring %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
fall %>%
  group_by(Site, Treatment) %>%
  filter(Season=='Spring') %>%
  count
# spring is missing a Lovelock control measurement, and contains one extra
# Lovelock full measuremen
spring %>%
  group_by(Site, Field, Treatment) %>%
  filter(Site=='Lovelock', Season=='Spring') %>%
  count
# one field 3 full measurement should be moved to field 3 control
spring %>%
  group_by(Field, Plot, Treatment) %>%
  filter(Site=='Lovelock', Field=='3', Season=='Spring') %>%
  count
# one field 3 plot 2 full measurement should be moved to field 3 plot 2 control
spring %>%
  filter(Site=='Lovelock', Field=='3', Season=='Spring') %>%
  mutate(VialID = paste(Vial, AllAph)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
# LLF3P2SHAM[?] is marked as full, but should probably be marked as control. The
# aphid count would not be out of place for a control measurement in this field

# Resolve error 8: move LLF3P2SHAM[?] -> Lovelock field 3, plot 2 control
spring %<>%
  mutate(Treatment = case_when(Vial=='LLF3P2SHAM[?]' & AllAph == 86 ~ 'Control',
                               Vial!='LLF3P2SHAM[?]' | AllAph != 86 ~
                                 as.character(Treatment)))

# Identify error 9
fall %>%
  group_by(Site, Field) %>%
  count
# Lovelock field 1 underweight, Lovelock includes one measurement with NA as
# field
fall %>%
  filter(Site=='Lovelock', Field==1 | is.na(Field)) %>%
  arrange(Plot)
# L7P3 full should probably be moved to Lovelock Field 1, plot 3 full
fall %>%
  filter(Site=='Lovelock', Field==1 | is.na(Field)) %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
# the L7P3 aphid count matches others from Lovelock field 1, plot 3
# Beware, this technically an empty string, not an NA value

# Resolve error 9: L7P3 full Lovelock field NA -> Lovelock field 1
fall %<>%
  mutate(Field = case_when(Field == '' ~ '1',
                           Field != '' ~ as.character(Field)))

# Identify error 10
fall %>%
  group_by(Site, Plot) %>%
  count
## Yerington plot 1 underweight, Yerington includes one measurement with NA as
## plot
fall %>%
  group_by(Field, Plot) %>%
  filter(Site=='Yerington') %>%
  arrange(Plot) %>%
  count
# Yerington field 2 plot 1 underweight, Yerington NA plot measurement is in
# field 1?
fall %>% filter(is.na(Plot))
fall %>% filter(Plot == '')
# Yikes - this is actually an empty string, not an NA value
# could YF1P7 possibly be Yerington Field 2, Plot 1?
fall %>% filter(Site=='Yerington', Field=='2', Plot=='1')
# Yerington field 2 plot 1 sham is missing, and the YF1P7 measurement is from a
# sham exclosure
# this suggests that YF1P7 should be moved to Yerington field 2 plot 1
fall %>% filter(Site=='Yerington', Field=='1' | Field=='2', Plot=='1' | is.na(Plot)) %>%
  mutate(VialID = paste(Vial, AllAph)) %>%
  ggplot(aes(x = VialID, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
# I can't say the data fit really makes sense - YF1P7 has an exceptionally high
# aphid count. I suppose it matches up with the expectation based on the high
# YF2P1 full count, but why is the YF2P1 (pre-) count so low? I'm going to go
# ahead and move YF1P7 into Yerington field 2 plot 1 anyway. This makes the most
# sense to me overall

# Resolve error 10: YF1P7 Yerington field 1 plot NA -> Yerington field 2 plot 1
fall %<>%
  mutate(Field = case_when(Plot == '' ~ '2',
                           !is.na(Plot) ~ as.character(Field)),
         Plot = case_when(Plot == '' ~ '1', !is.na(Plot) ~ as.character(Plot)))

# Identify error 11
spring %>%
  group_by(Field, Plot, Treatment, Season) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
# field 3 plot 1 pre- underweight, field 3 plot 3 pre- overweight
spring %>%
  group_by(Site, Field, Plot, Treatment, Season) %>%
  filter(Season=='Spring', Field=='3', Plot=='1' | Plot=='3') %>%
  arrange(Field, Plot, Treatment) %>%
  count
# the problem is contained to Yerington
spring %>%
  group_by(Field, Plot, Season) %>%
  filter(Season=='Spring',
         Site=='Yerington',
         Field=='3',
         Plot=='1' | Plot=='3') %>%
  arrange(Field, Plot, Treatment)
# could be that YF2.2(Location3) is the one that needs to be moved to plot 1.
# this makes sense that YF2.2(Location3)P3 specifically has P3 noted as the
# plot, while YF2.2(Location3) lacks any specific plot designation it seems
# likely that I would forget to record the plot number on the first plot, but
# not on the third plot
spring %>%
  group_by(Field, Plot, Season) %>%
  filter(Season=='Spring',
         Site=='Yerington',
         Field=='3',
         Plot=='1' | Plot=='3') %>%
  ggplot(aes(x = Vial, y = AllAph, fill = Treatment)) +
  geom_bar(stat = 'identity')
# the data here is variable - in some treatments, plot 1 is higher; in others,
# plot 3 is higher this doesn't inform which of the field 3 plot 3 measurements
# should be moved aphid counts are generally very low anyway, so some
# stochasticity is expected all I have to go on is the way the labels were
# recorded, but that is already enough to convince me

# Resolve error 11: YF2.2(Location3) plot 3 -> plot 1
spring %<>%
  mutate(Plot = case_when(Vial=='YF2.2(Location3)' ~ '1',
                          Vial!='YF2.2(Location3)' ~ as.character(Plot)))

# final checks
spring %>%
  group_by(Field, Plot, Treatment, Season) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
fall %>%
  group_by(Field, Plot, Treatment, Season) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
# spring data looks complete
fall %>%
  group_by(Field, Plot, Treatment, Season) %>%
  filter(Season=='Fall') %>%
  arrange(Field, Plot) %>%
  count
spring %>%
  group_by(Field, Plot, Treatment, Season) %>%
  filter(Season=='Fall') %>%
  arrange(Field, Plot) %>%
  count
# fall data looks complete

# merge spring and fall data into a single tibble
data <- rbind(spring %>% select(-Notes),
              fall %>% select(-Number)) %>%
  mutate(Season = as_factor(Season))
data %>%
  group_by(Field) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
# field 4 should move into field 2 across the entire dataset
data %<>%
  mutate(Field = case_when(Field=='4' ~ '2', Field!='4' ~ as.character(Field)))
data %>%
  group_by(Site, Field, Plot) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
# balanced!
data %>%
  group_by(Site, Field, Treatment) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot) %>%
  count
# balanced!
data %>%
  group_by(Field, Plot, Treatment) %>%
  filter(Season=='Spring') %>%
  arrange(Field, Plot, Treatment) %>%
  count
# remove 'spring' and 'fall' tibbles
rm(spring, fall)

## check veg data ####
# Identify error 1
str(vegdata)
levels(vegdata$species)
vegdata %>%
  mutate(species = fct_relevel(species, sort)) %>%
  group_by(species) %>%
  count
# many species are incorrectly ID'ed or unknown
# joining the updated species id table can help
new_ids %>% group_by(orig_code) %>% count
vegdata %>%
  left_join(new_ids,  by = c("species" = "orig_code")) %>%
  mutate(species = case_when(is.na(new_code) ~ paste(species),
                             !is.na(new_code) ~ paste(new_code))) %>%
  mutate(species = case_when(species == '' ~ paste(NA), # also clean up NAs
                             species == 'NA' ~ paste(NA),
                             is.na(species) ~ paste(NA),
                             species == 'none' ~ paste(NA),
                             !is.na(species) ~ paste(species))) %>%
  mutate(species = fct_relevel(species, sort)) %>%
  group_by(species) %>%
  count
# correct error 1
# join updated_species_id.csv
vegdata %<>%
  left_join(new_ids,  by = c("species" = "orig_code")) %>%
  mutate(species = case_when(is.na(new_code) ~ paste(species),
                             !is.na(new_code) ~ paste(new_code))) %>%
  mutate(species = case_when(species == '' ~ NA_character_, # also clean up NAs
                             species == 'NA' ~ NA_character_,
                             is.na(species) ~ NA_character_,
                             species == 'none' ~ NA_character_,
                             !is.na(species) ~ paste(species)))
# rm 'new_ids'
rm(new_ids)

# identify error 2
vegdata %>%
  group_by(species) %>%
  summarize(mean_cover = mean(cover), .groups = 'keep') %>%
  arrange(desc(mean_cover)) %>%
  filter(species %in% c('MESA', 'ALFA'))
# some MESA recorded as ALFA
# correct error 2
vegdata %<>%
  mutate(species = case_when(species == 'ALFA' ~ paste('MESA'),
                             species != 'ALFA' ~ paste(species)))
# identify error 3
vegdata %>%
  group_by(site) %>%
  summarize(count = n_distinct(waypoint))
# note that we do not expect sites to have even numbers of waypoints (some were
# inaccessible)
vegdata %>%
  filter(is.na(count), !is.na(species))
# missing count values are also expected
vegdata %>%
  filter(is.na(cover), !is.na(species))

# missing cover values are NOT expected
# note that ALCE, MESA are ag plots that were observed from a distance. Cover is
# likely high. In othen instances, cover was likely negligible.
# correct error 2
# build table of replacement cover values
new_cover <- vegdata %>%
  filter(is.na(cover), !is.na(species)) %>%
  mutate(newcover = c(1, 50, 50, 50, 50, 1, 1, 50, 50, 1, 1, 1)) %>%
  select(survey_date, site, waypoint, species, newcover)
# join new cover
vegdata %<>%
  left_join(new_cover, by = c('survey_date', 'site', 'waypoint', 'species')) %>%
  mutate(cover = case_when(is.na(cover) ~ newcover, !is.na(cover) ~ cover)) %>%
  # drop unwanted cols
  select(survey_date, observers, site, waypoint, species, count, cover)

# remove new_cover var
rm(new_cover)

# identify error 5
vegdata %>%
  filter(site == 'Lovelock', waypoint == 26) %>%
  arrange(species)
# duplicate BASC5 counts
# correct error 5
# remove the offending row
vegdata %<>% filter(!(site == 'Lovelock' &
                        waypoint == 26 &
                        species == 'BASC5' &
                        count == 3))

# identify error 6
vegdata %>%
  filter(site == 'Fallon', waypoint == 7) %>%
  arrange(species)
# duplicate CHVI4 counts
# correct error 6
# remove the offending row
vegdata %<>% filter(!(site == 'Fallon' &
                        waypoint == 7 &
                        species == 'CHVI4' &
                        count == 1))

# identify error 7
vegdata %>%
  filter(site == 'Yerington', waypoint == 7) %>%
  arrange(species)
# duplicate HAGL counts
# correct error 7
# remove the offending row
vegdata %<>% filter(!(site == 'Yerington' &
                        waypoint == 7 &
                        species == 'HAGL' &
                        count == 3))

## check location data ####
str(rawlocs)
# col spec looks good
rawlocs %>% group_by(site, CID) %>% count
# uneven weights are expected
# everything looks ok

# check veg plots cover data
veg_survey_cover %>% group_by(site) %>% count
# uneven weights are expected
# not really much else to check

## check raw cover data ####
# checked after wrangling

# data wrangling ####
## insect counts ####
# lengthen
data_long <- pivot_longer(data,
                          cols = Arachnida:NonAcy,
                          names_to = 'Taxa',
                          values_to = 'Density') %>%
  mutate(Taxa = as_factor(Taxa))
# split off and re-join 'Pre' data; select relevant columns; widen
data <- right_join(data_long %>% filter(Treatment != 'Pre-'),
                   data_long %>% filter(Treatment == 'Pre-'),
                   by = c('Site','Field','Plot','Taxa','Season')) %>%
  # rename cols where necessary
  select(Site:Plot,
         Treatment = Treatment.x,
         Taxa,
         Density = Density.x,
         Pre_Density = Density.y,
         Season) %>%
  # widen to put each taxon in a unique col
  pivot_wider(names_from = Taxa, values_from = Density:Pre_Density) %>%
  # relevel factors to reflect natural order of treatments and seasons
  mutate(Treatment = fct_relevel(Treatment, c('Control', 'Sham', 'Full')),
         Season = fct_relevel(Season, c('Spring','Fall'))) %>%
  # create 'Non_Acy' cols
  mutate(Pre_Density_Non_Acy = Pre_Density_Aphis + Pre_Density_Therioaphis,
         Density_Non_Acy = Density_Aphis + Density_Therioaphis)
# final check
str(data)

## locations #####
# wrangle field location data
locs <- rawlocs %>%
  filter(CID != 7) %>% # remove non-survey points
  mutate(id = paste0(site,
                     gsub('^\\D*(\\d*)(ALT)?\\D*$',
                          '\\1\\2',
                          .$Name))) %>% # parse site col
  select(id, field_id, long = POINT_X, lat = POINT_Y) %>% # drop unused cols
  left_join(veg_survey_classes %>%
              mutate(id = paste0(site, plotnum)), # join veg_survey_classes from ee data
            by = 'id') %>%
  select(-name) %>%  # drop 'name' col
  relocate(id, site, plotnum, type, lat, long, field_id) %>%  # reorder cols
  left_join(veg_survey_cover %>% # tidy up veg_survey_cover before joining
              unite(id, site, plotnum, sep = '') %>%
              select(id, K3LSfall:XsentSummer), by = 'id')

# rm 'rawlocs,' 'veg_survey_classes'
rm(rawlocs, veg_survey_classes, veg_survey_cover)

## veg data #####
# make season and site ID column
vegdata %<>%
  # use getSeason() to assign 'season' based on 'survey_date'
  mutate(season = getSeason(survey_date)) %>%
  unite(id, site, waypoint, sep = '', remove = FALSE) %>% # create id col
  left_join(locs, by = 'id') %>% # join location data
  relocate(id, lat, long, site = site.x, waypoint, type, field_id, season,
           survey_date, observers, species, count, cover) # reorder cols

# two plots have no location data
# probably recording errors that aren't resolvable
# Yerington 8B is a duplicate, throw out
# Fallon 9 is probably a recording error that isn't resolvable, throw out
vegdata %<>% filter(!is.na(lat))

## landcover data ####
# build list of col names
klasses <- map(seq(0, 11), toString)
distanceWeight <- c('no',
                    'const',
                    'sig1',
                    'sig2',
                    'sig3',
                    'sig4',
                    'sig5')
testnames <- matrix(NA, 12, 7)
for (i in 1:length(distanceWeight)) {
  for (j in 1:length(klasses))
  testnames[[j, i]] <- paste0(klasses[[j]], distanceWeight[[i]])
}
testnames
# build list of site names
site <- c('Minden',
              'Minden',
              'Minden',
              'Lovelock',
              'Lovelock',
              'Lovelock',
              'Fallon',
              'Fallon',
              'Fallon',
              'Yerington',
              'Yerington',
              'Yerington')
field <- c('01',
               '02',
               '03',
               '03', # note non-numeric order here. this is correct.
               '01',
               '02',
               '01',
               '02',
               '03',
               '01',
               '02',
               '03')
sitename = c()
for (i in 1:length(field)) {
  sitename[[i]] = paste0(site[[i]], field[[i]])
}
sitename
# tidy
tidyLandcover <- landcover %>%
  # split cols
  separate(col = dataSeries,
           into = testnames,
           sep = ",") %>%
  # add siteId
  mutate(siteId = unlist(sitename), .before = everything()) %>%
  # begin to transpose
  pivot_longer(-siteId, 'var') %>%
  # parse cells
  separate(col = 'value',
           into = c('class', 'value'),
           sep = '=') %>%
  rowwise %>%
  mutate(class = paste0('class', str_extract(class, '[:digit:]+')),
         distanceWeight = str_extract(var, '[:alpha:]+[:digit:]?'),
         value = parse_number(map_chr(str_extract_all(value,
                                                      '[:digit:]+|\\.|-|E'),
                                      ~ str_c(.x, collapse = '')))) %>%
  # finish transpose
  pivot_wider(c(siteId, distanceWeight), names_from = class) %>%
  # make site and field id
  separate(col = siteId,
           into = c('site', 'field'),
           sep = -2,
           remove = FALSE)
View(tidyLandcover)
summary(tidyLandcover)

# completeness/accuracy check: check areaScore * class (should be equal within
# image sets and decay functions)

# get sum of area scores for all classes
df <- tidyLandcover %>%
  mutate(sumScore = sum(across(where(is.double)))) %>%
  group_by(site, distanceWeight) %>%
  summarize(totalSum = sum(sumScore), .groups = 'keep')
# filter to a specific decay function
no <- df %>%
  filter(distanceWeight == 'no') %>%
  arrange(desc(totalSum))
# divide largest area by smallest area
noRatio <- head(no$totalSum, n = 1)/tail(no$totalSum, n = 1)

const <- df %>%
  filter(distanceWeight == 'const') %>%
  arrange(desc(totalSum))

constRatio <- head(const$totalSum, n = 1)/tail(const$totalSum, n = 1)

sig1 <- df %>%
  filter(distanceWeight == 'sig1') %>%
  arrange(desc(totalSum))

sig1Ratio <- head(sig1$totalSum, n = 1)/tail(sig1$totalSum, n = 1)

sig2 <- df %>%
  filter(distanceWeight == 'sig2') %>%
  arrange(desc(totalSum))

sig2Ratio <- head(sig2$totalSum, n = 1)/tail(sig2$totalSum, n = 1)

sig3 <- df %>%
  filter(distanceWeight == 'sig3') %>%
  arrange(desc(totalSum))

sig3Ratio <- head(sig3$totalSum, n = 1)/tail(sig3$totalSum, n = 1)

sig4 <- df %>%
  filter(distanceWeight == 'sig4') %>%
  arrange(desc(totalSum))

sig4Ratio <- head(sig4$totalSum, n = 1)/tail(sig4$totalSum, n = 1)

sig5 <- df %>%
  filter(distanceWeight == 'sig5') %>%
  arrange(desc(totalSum))

sig5Ratio <- head(sig5$totalSum, n = 1)/tail(sig5$totalSum, n = 1)

out <- c(noRatio, constRatio, sig1Ratio, sig2Ratio,
         sig3Ratio, sig4Ratio, sig5Ratio)

out # perfection!!!

# remove misc error-checking vars from env
rm(df, no, const, sig1, sig2, sig3, sig4, sig5,
   noRatio, constRatio, sig1Ratio, sig2Ratio, sig3Ratio, sig4Ratio, sig5Ratio,
   distanceWeight, out, field, site)

## calculate vegdata diversity metrics ####
# these metrics are all based on cover, since this is the more complete data
### make subplot-level diversity table ####
subplot_divdata <- vegdata %>%
  mutate(id_season = paste0(id, season)) %>%
  select(id_season, id, season, species, cover) %>%
  pivot_wider(names_from = species, values_from = cover, values_fill = 0) %>%
  select(-'NA')
# calculate subplot-level shannon diversity
shan <- subplot_divdata %>%
  select(-id_season, -id, -season) %>%
  diversity(index = 'shannon')
# calculate subplot-level simpson diversity
simp <- subplot_divdata %>%
  select(-id_season, -id, -season) %>%
  diversity(index = 'simpson')
# calculate subplot-level species richness
rich <-subplot_divdata %>%
  select(-id_season, -id, -season) %>%
  specnumber()
# cbind diversity metrics and subplot info
subplot_divdata %<>%
  ungroup() %>%
  mutate(shan = all_of(shan), simp = all_of(simp), rich = all_of(rich))

### join diversity table and location data ####
vegPlots <- subplot_divdata %>%
  left_join(locs, by = 'id') %>%
  select(-id) %>% # drop simple (no season) id since this is no longer needed
  mutate(type = as_factor(type)) %>% # covert 'type' to factor
  # reorder all cols, fix capitalization
  relocate(id = id_season, field_id, site, plotnum, type, season, lat, long, shan, simp, rich)

# rm helper vars
rm(subplot_divdata, shan, simp, rich, locs)

### make site * cover class level diversity table ####
# this is based on means, which is arbitrary (could also be median, etc.) I'm
# not worried about it right now because I'm not ever sure how this data will be
# used
vegSites <- vegdata %>%
  unite(id, id, season, remove = FALSE) %>%
  select(site, id, species, cover) %>%
  pivot_wider(names_from = species, values_from = cover, values_fill = 0) %>%
  select(-'NA') %>%
  group_by(site) %>%
  summarise(across(ACHY:ZEMA, ~ mean(.x))) # using mean()
# calculate simpson diversity
shan <- vegSites %>%
  ungroup() %>%
  select(-site) %>%
  diversity(index = 'shannon')
# calculate simpson diversity
simp <- vegSites %>%
  ungroup() %>%
  select(-site) %>%
  diversity(index = 'simpson')
# calculate species richness
rich <- vegSites %>%
  ungroup() %>%
  select(-site) %>%
  specnumber
# cbind diversity metrics and subplot info
vegSites %<>%
  ungroup() %>%
  mutate(shan = all_of(shan), simp = all_of(simp), rich = all_of(rich)) %>%
  relocate(site, shan, simp, rich)

# ee data not yet added to vegSites
# not sure this would be useful anyway
# rm helper vars
rm(shan, simp, rich, vegdata)

# define convenience identifiers ####
# move this to main script
# define lists of predator and aphid taxa
predlist <- c('Arachnida','Coccinellidae','Ichneumonidae',
              'Nabis', 'Geocoris', 'Anthocoridae')
aphlist <- c('Acyrthosiphon', 'Aphis', 'Therioaphis', 'AllAph', 'NonAcy')

# define predator- and aphid-only tibbles
# filter data to include only predators
preddata <- right_join(data_long %>% filter(Treatment != 'Pre-'),
                       data_long %>% filter(Treatment == 'Pre-'),
                       by = c('Site','Field','Plot','Taxa','Season')) %>%
  select(Site:Plot,
         Treatment = Treatment.x,
         Taxa,
         Density = Density.x,
         Pre_Density = Density.y,
         Season) %>%
  filter(Taxa %in% predlist)

aphdata <- right_join(data_long %>% filter(Treatment != 'Pre-'),
                      data_long %>% filter(Treatment == 'Pre-'),
                      by = c('Site','Field','Plot','Taxa','Season')) %>%
  select(Site:Plot,
         Treatment = Treatment.x,
         Taxa,
         Density = Density.x,
         Pre_Density = Density.y,
         Season) %>%
  filter(Taxa %in% aphlist)

# filter out 'Full' data ####
# remove 'Full' treatments from all data
data %<>% filter(Treatment != 'Full')
data_long %<>% filter(Treatment != 'Full')

# preview tidy data ####
# print a summary of each final dataset
data # main table of insect count data
data_long # insect data in long format
tidyLandcover # landcover areaScores for each field
vegPlots # final plot-level veg data from surveys
vegSites # final site-level veg data from surveys

# export tidy data ####
# build list of data to export
tidy_data <- list(data, data_long, tidyLandcover, vegPlots, vegSites)
data_names <- list('data', 'data_long', 'landcover', 'vegPlots', 'vegSites')
# export
walk2(tidy_data, data_names, ~write_csv(.x, paste0('tidy_data/', .y, ".csv")))

