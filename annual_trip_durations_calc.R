library(magrittr); library(dplyr)

load("Data_outputs/all_gps_data.RData")
load("Data_inputs/ind_meta.RData")

all_ind.meta %<>% mutate(species = tolower(species))

all_gps.df <- merge(all_gps.df, all_ind.meta, by= c("ring", "species"))

annual_trip.durs <- all_gps.df %>% 
  group_by(ring, boutID) %>% slice(1) %>% 
  ungroup() %>%
  group_by(colony, species, season, sex, phase) %>% 
  summarise(mean_trip.hrs = mean(duration.hours, na.rm = T),
            med_trip.hrs = median(duration.hours, na.rm =T),
            sd_trip.hrs = sd(duration.hours, na.rm = T),
            n_birds = n_distinct(ring))

write.csv(annual_trip.durs, file = "Data_outputs/annual_trip_durations.csv", row.names = F)
