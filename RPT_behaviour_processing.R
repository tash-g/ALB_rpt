## ---------------------------
##
## Script name: individual_repeatability_processing
##
## Purpose of script: Process GLS & GPS data for individual repeatability estimates.
## NOTE This script interfaces with the main ALB project.
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-10-14
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
#source("ALB_FOR_functions.R")
#source("ALB_data_processing_functions.R")
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "readxl", "sf", "ggplot2", "data.table",
              "momentuHMM")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], dependencies =T)
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)
select <- dplyr::select



# Get sample sizes --------------------------------------------------------

### GLS ----------
# Change directory to access GLS data
setwd("C:/Users/gilli/OneDrive - The University of Liverpool/Liverpool postdoc/ALB_foraging/ALB_foraging_proj")

spec_col <- c("bba_birdis", "bba_kerguelen", "waal_birdis", "waal_crozet")
sample_gls.list <- list()
runs_gls.list <- list()

for (i in 1:length(spec_col)) {
  
  load(paste0("Data_outputs/", spec_col[i], "_gls_labelled_subset.RData"))
  
  seasons.table <- gls_labelled.df %>% 
    group_by(ring) %>%
    summarise(n_seasons = n_distinct(season)) %>%
    mutate(spec_col = spec_col[i]) %>%
    separate(spec_col, into = c("species", "colony"), sep = "_")
  
  sample_gls.list[[i]] <- seasons.table
  
  
  #### Duration per individual --------
  
  # Get simplified daily dataframe
  gls_daily <- gls_labelled.df %>% 
    mutate(date = as.Date(datetime)) %>%
    group_by(ring, season, date) %>%
    arrange(date) %>%
    select(ring, season, date) %>%
    distinct()
  
  # Find number of continuous runs for each individual
  gls_daily %<>%
    arrange(ring, season, date) %>%
    group_by(ring, season) %>%
    mutate(date_num = as.integer(date),
           gap = date_num - lag(date_num, default = first(date_num)),
           new_run = gap != 1,
           # Create run ID by cumulative sum of new runs (start a new run when gap != 1)
           run_id = cumsum(new_run) + 1 )
  
  # Summarise to get run length and start/end dates per run per bird-colony
  runs_summary <- gls_daily %>%
    group_by(ring, season, run_id) %>%
    summarise(
      run_start = min(date),
      run_end = max(date),
      run_length = n() ) %>%
    select(-run_id) %>%
    mutate(spec_col = spec_col[i]) %>%
    separate(spec_col, into = c("species", "colony"), sep = "_")
 
  runs_gls.list[[i]] <- runs_summary
  
}

sample_sizes.gls <- do.call("rbind", sample_gls.list)

sample_sizes.gls %>%
  mutate(n_seasons_group = case_when(
      n_seasons == 1 ~ "1_season",
      n_seasons == 2 ~ "2_seasons",
      n_seasons == 3 ~ "3_seasons",
      n_seasons >= 4 ~ "4+_seasons",
      TRUE ~ NA_character_) ) %>%
  group_by(species, colony, n_seasons_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = n_seasons_group,
    values_from = n,
    values_fill = 0) %>%
  mutate( total = `1_season` + `2_seasons` + `3_seasons` + `4+_seasons`)


runs_gls <- do.call("rbind", runs_gls.list)

runs_gls %>%
  group_by(ring) %>%
  mutate(n_total = sum(run_length)) %>%
  group_by(species, colony) %>%
  summarise(mean_dur = mean(n_total),
            med_dur = median(n_total),
            sd_dur = sd(n_total),
            n_birds = n_distinct(ring))

### GPS ---------

setwd("C:/Users/gilli/OneDrive - The University of Liverpool/Liverpool postdoc/ALB_foraging/ALB_rpt")

sample_gps.list <- list()
runs_gps.list <- list()

for (i in 1:length(spec_col)) {
  
  load(paste0("Data_inputs/", spec_col[i], "_gps_labelled_subset.RData"))
  
  seasons.table <- gps_labelled.df %>% 
    group_by(ring) %>%
    summarise(n_seasons = n_distinct(season)) %>%
    mutate(spec_col = spec_col[i]) %>%
    separate(spec_col, into = c("species", "colony"), sep = "_")
  
  sample_gps.list[[i]] <- seasons.table
  
  #### Duration per individual --------
  
  gps_daily <- gps_labelled.df %>% 
    mutate(date = as.Date(datetime)) %>%
    group_by(ring, season, date) %>%
    arrange(date) %>%
    select(ring, season, date) %>%
    distinct() 
  
  # Find number of continuous runs for each individual
  gps_daily %<>%
    arrange(ring, season, date) %>%
    group_by(ring, season) %>%
    mutate(date_num = as.integer(date),
           gap = date_num - lag(date_num, default = first(date_num)),
           new_run = gap != 1,
           # Create run ID by cumulative sum of new runs (start a new run when gap != 1)
           run_id = cumsum(new_run) + 1 )
  
  # Summarise to get run length and start/end dates per run per bird-colony
  runs_summary <- gps_daily %>%
    group_by(ring, season, run_id) %>%
    summarise(
      run_start = min(date),
      run_end = max(date),
      run_length = n() ) %>%
    select(-run_id) %>%
    mutate(spec_col = spec_col[i]) %>%
    separate(spec_col, into = c("species", "colony"), sep = "_")
  
  runs_gps.list[[i]] <- runs_summary
  
}

sample_sizes.gps <- do.call("rbind", sample_gps.list)

sample_sizes.gps %>%
  mutate(n_seasons_group = case_when(
    n_seasons == 1 ~ "1_season",
    n_seasons == 2 ~ "2_seasons",
    n_seasons == 3 ~ "3_seasons",
    n_seasons >= 4 ~ "4+_seasons",
    TRUE ~ NA_character_) ) %>%
  group_by(species, colony, n_seasons_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = n_seasons_group,
    values_from = n,
    values_fill = 0) %>%
  mutate( total = `1_season` + `2_seasons` + `3_seasons` + `4+_seasons`)


runs_gps <- do.call("rbind", runs_gps.list)

runs_gps %>%
  group_by(ring) %>%
  mutate(n_total = sum(run_length)) %>%
  group_by(species, colony) %>%
  summarise(mean_dur = mean(n_total),
            med_dur = median(n_total),
            sd_dur = sd(n_total),
            n_birds = n_distinct(ring))

# ______________________________ ####
# ~ DATA PROCESSING ~ ###########################################################

# +++++++++++++++++++++++++++++ ####

# CUT GLS DATA TO BREEDING ---------------------------------------------------------

# Change directory to access GLS data
setwd("C:/Users/gilli/OneDrive - The University of Liverpool/Liverpool postdoc/ALB_foraging/ALB_foraging_proj")

spec_col <- c("bba_birdis", "bba_kerguelen", "waal_birdis", "waal_crozet")

for (i in 1:length(spec_col)) {
  
  load(paste0("Data_inputs/", spec_col[i], "_gls_labelled.RData"))
  
  gls_labelled.df %<>%
    mutate(ring = toupper(ring),
           ring = gsub(" ", "", ring),
           ringYr = paste(ring, season, sep = "_")) %>%
    filter(!is.na(ring))
  
  
  # Load meta data -------------------------------------------------------------
  
  meta_source <- ifelse(grepl("birdis", spec_col[i]), "BAS", "Chize")
  all_demo <- loadRData(paste0("Data_inputs/", meta_source, "_demo_complete.RData"))
  
  all_demo %<>% mutate(species = tolower(species)) %>% filter(species == my_species) 
  
  ## Make dataframe long
  demoA <- all_demo %>% select(-c(bird2_sex, bird2_birthYr)) %>% 
    dplyr::rename(ring = bird1, sex = bird1_sex, birthYr = bird1_birthYr, partner = bird2) %>%
    filter(!is.na(ring)) 
  
  demoB <- all_demo %>% select(-c(bird1_sex, bird1_birthYr)) %>% 
    dplyr::rename(ring = bird2, sex = bird2_sex, birthYr = bird2_birthYr, partner = bird1) %>%
    filter(!is.na(ring)) 
  
  my_demo <- rbind(demoA, demoB)
  my_demo %<>% group_by(ring, season) %>% tidyr::fill(everything(), .direction = 'updown') %>% distinct()
  
  
  
  # Manual breeding dates ------------------------------------------------------
  
  if (spec_col[i] == "waal_birdis") {
    
    median_lay <- "12-01" 
    median_hatch <- "03-10"
    
  } else if (spec_col[i] == "bba_birdis") {
    
    median_lay <- "10-13"
    median_hatch <- "01-05"
    
  } else if (spec_col[i] == "waal_cro") {
    
    median_lay <- "12-01"
    median_hatch <- "03-22"
    
  } else if (spec_col[i] == "bba_kerguelen") {
    
    median_lay <- "09-25"
    median_hatch <- "12-15"
  }
  
  
  # Match in metadata and filter to breeding -------------------------------------
  
  setDT(my_demo)
  my_demo %<>% select(ring, season, rs, lay_date, hatch_date, fail_date) %>%
    mutate(season = as.numeric(season)) %>% distinct()
  my_demo <- my_demo[, .SD[sample(.N, 1)], by = .(ring, season)]
  
  setDT(gls_labelled.df)
  
  # Merge using data.table
  gls_labelled.df <- merge(gls_labelled.df, my_demo, by = c("ring", "season"), all.x = T)
  
  ## Cut to breeding
  hatch_time = ifelse(my_species == "bba", 101, 121)
  brood_time = ifelse(my_species == "bba", 38, 48) # brooding ends 22/01
  
  gls_labelled.df %<>% rename(boutID = tripID)
  
  gls_labelled.df <- cut_to_breeding(gls_labelled.df, median_lay, median_hatch, hatch_time, brood_time)
  
  save(gls_labelled.df, file = paste0("Data_outputs/", spec_col[i], "_gls_labelled_subset.RData"))

}
  
# +++++++++++++++++++++++++++++ ####

# PROCESS GPS DATA -------------------------------------------------------------

load(paste0("Data_inputs/", my_species, "_", colony_exp, "_gps_labelled_subset.RData"))

## Calculate landings using speed filter ---------------------------------------

# Remove trips with unrealistic speed
gps_labelled.df %<>% filter(calc_speed < 26)

# Deal with duplicates
gps_labelled.df %<>%
  group_by(ID) %>%
  mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>%
  distinct() %>%
  arrange(datetime)

# Label landings & transit
gps_labelled.df %<>% 
  mutate(behav_state = ifelse(calc_speed < 2.7, "rest", "aloft"),
         prev_state = lag(behav_state),
         landing = ifelse(behav_state == "aloft" & prev_state == "rest",
                          "take-off",
                          ifelse(behav_state == "rest" & prev_state == "aloft",
                                 "landing", "none")),
         landing = ifelse(is.na(landing), "none", landing),
         behav_state = ifelse(calc_speed > 8.3, "transit", behav_state))


# Calculate distance to next point
gps_labelled.df %<>%
  group_by(ring, boutID) %>%
  mutate(dist_next = gcd.hf(longitude, latitude, lead(longitude), lead(latitude))) %>%
  ungroup()

## Summarise behavioural metrics ----------------------------------------------

# gps_summary <- gps_labelled.df %>%
#   group_by(ring, season, boutID) %>%
#   filter(loc == "trip") %>%
#   summarise(ring = ring[1],
#             season = season[1],
#             ringYr = ringYr[1],
#             start_date = min(datetime),
#             end_date = max(datetime),
#             duration.hours = as.numeric(difftime(end_date, start_date, units = "hours") + 0.17),
#             duration.days = duration.hours / 24,
#             n_landings = sum(landing != "none"),
#             transit_time.hrs = (sum(behav_state == "transit") * 10)/60,
#             rest_time.hrs = (sum(behav_state == "rest") * 10)/60,
#             transit_time.pct = (transit_time.hrs / duration.hours) * 100,
#             rest_time.pct = (rest_time.hrs / duration.hours) * 100,
#             total_distance.km = sum(dist_next),
#             max_distance.km = max(distances) / 1000 )

## Full dataset to try integrating with GLS
gps_full <- gps_labelled.df %>%
  group_by(ring, season, boutID) %>%
  mutate(start_date = min(datetime),
         end_date = max(datetime),
         duration.hours = as.numeric(difftime(end_date, start_date, units = "hours") + 0.17),
         duration.days = duration.hours / 24,
         n_landings = sum(landing != "none"),
         transit_time.hrs = (sum(behav_state == "transit") * 10)/60,
         rest_time.hrs = (sum(behav_state == "rest") * 10)/60,
         transit_time.pct = (transit_time.hrs / duration.hours) * 100,
         rest_time.pct = (rest_time.hrs / duration.hours) * 100,
         total_distance.km = sum(dist_next),
         max_distance.km = max(dist_col) / 1000 ) %>%
  select(-c(start_date, end_date))



# +++++++++++++++++++++++++++++ ####

# PROCESS GLS DATA WITH GPS DATA ===============================================

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_gls_labelled_subset.RData"))

all_birds <- unique(c(unique(gls_labelled.df$ringYr),
                      unique(gps_full$ringYr)))
all_trips.list <- vector(mode = "list", length = length(all_birds))

pb <- txtProgressBar(min = 0, max = length(all_birds), style = 3)

for (i in 1:length(all_birds)) {
  setTxtProgressBar(pb, i)
  
  mygps <- subset(gps_full, ringYr == all_birds[i]) %>%
    mutate(immersion = NA, dur.mins = NA, 
           dryWet = ifelse(behav_state == "rest", "wet", "dry"),
           source = "GPS") %>%
    ungroup() %>%
    select(c(ring, season, longitude, latitude, datetime, ringYr, phase, loc, rs, 
             lay_date, hatch_date, hatch_code, fail_date, 
             n_landings, behav_state, total_distance.km,
             max_distance.km, immersion, dryWet, source)) 
  
  mygls <- subset(gls_labelled.df, ringYr == all_birds[i]) %>%
    mutate(longitude = NA, latitude = NA, n_landings = NA, total_distance.km = NA,
           max_distance.km = NA, source = "GLS",
           behav_state = ifelse(immersion <= 3, "transit", 
                                ifelse(immersion >= 197, "rest", "aloft")),
           datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
    select(c(ring, season, longitude, latitude, datetime, ringYr, phase, loc, rs, 
             lay_date, hatch_date, hatch_code, fail_date, n_landings, behav_state, 
             total_distance.km, max_distance.km, immersion, dryWet, source)) %>%
    relocate(ring, season, longitude, latitude, datetime, ringYr, phase, loc, rs, 
             lay_date, hatch_date, hatch_code, fail_date, n_landings, 
             behav_state, total_distance.km, max_distance.km, immersion, dryWet, source)

  if(nrow(mygps) > 0 & nrow(mygls) > 0) {
    # Round to nearest minute and where duplicates, retain GPS
    combined <- rbind(mygps, mygls) %>% 
      mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
      arrange(datetime) %>%
      mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>%
      group_by(datetime) %>%  
      filter(source == "GPS" | n() == 1) %>%  
      ungroup() %>% 
      mutate(agreement = ifelse(loc == lead(loc), "Y", "N"),
             loc2 = ifelse(agreement == "N" & source == "GLS", lead(loc), loc),
             loc2 = ifelse(is.na(loc2), loc, loc2)) %>%
      select(-c(loc, agreement, source)) %>%
      rename(loc = loc2)
    
  } else if (nrow(mygps) > 0 & nrow(mygls) == 0) { 
    combined <- mygps %>%
      arrange(datetime) %>%
      mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>% 
      distinct() } else { combined <- mygls %>%
    arrange(datetime) %>%
    mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>% 
    group_by(datetime) %>%
    slice_max(immersion, with_ties = FALSE) %>%
    ungroup() %>%
    distinct() }
  
  ### Extract trips ------------------------------------------------------------
  
  combined %<>% arrange(datetime)
  combined$tripID = rep(1:length(rle(combined$loc)$values), rle(combined$loc)$lengths)
  
  my_trips <- combined %>% 
    group_by(ring, season, tripID) %>%
    mutate(start = datetime[1],
           end = datetime[n()],
           duration.mins = as.numeric(difftime(end, start, unit = "mins"))) %>%
    arrange(start) %>%
    filter(!is.na(datetime) & !is.na(end))
  
  ## Get rid of very short trips/col visits
  if (all(my_trips$duration.mins < 480)) next
  
  max_tries <- 5
  tries <- 0
  
  while(any(my_trips$duration.mins <= 480) & tries <= max_tries) {
    
    tries = tries + 1
    
    my_trips$loc <- ifelse(my_trips$duration.mins <= 480, 
                           ifelse(my_trips$loc == "trip", "col", "trip"),
                           my_trips$loc)
    my_trips$tripID <- rep(1:length(rle(my_trips$loc)$values), rle(my_trips$loc)$length)
    
    my_trips %<>%
      group_by(ring, season, tripID) %>%
      mutate(start = datetime[1],
             end = datetime[n()],
             duration.mins = as.numeric(difftime(end, start, unit = "mins"))) %>%
      arrange(start)
    
  }
  
  # Summarise behaviours
  my_trips %<>%
    group_by(tripID) %>%
    filter(loc == "trip") %>%
    mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
           prev_state = lag(dryWet),
           bout_time = as.numeric(difftime(datetime, lag(datetime), units = "mins"))) %>%
    dplyr::summarise(species = my_species,
                     colony = colony_exp,
                     ring = ring[1],
                     season = season[1],
                     phase = phase[1],
                     start_date = datetime[1],
                     end_date = max(datetime, na.rm = T),
                     loc = loc[1],
                     # Duration
                     duration.hrs = as.numeric(difftime(end_date, start_date, units = "hours") + 0.17),
                     duration.days = duration.hrs / 24,
                     duration.hrs.var = var(duration.hrs),
                     # Transitions
                     transitions = sum(dryWet != prev_state, na.rm = TRUE),  # Count state transitions
                     transitions.var = var(transitions),
                     landing = sum(dryWet == "dry" & prev_state == "wet", na.rm = TRUE),
                     takeoff = sum(dryWet == "wet" & prev_state == "dry", na.rm = TRUE),
                     # Behaviours
                     transit_time.hrs = sum(bout_time[behav_state == "transit"], na.rm = T) / 60,
                     rest_time.hrs = sum(bout_time[behav_state == "rest"], na.rm = T) / 60,
                     transit_time.pct = transit_time.hrs/duration.hrs * 100,
                     rest_time.pct = rest_time.hrs/duration.hrs * 100,
                     transit_to_rest.ratio = transit_time.hrs/rest_time.hrs,
                     # Distances
                     total_distance.km = sum(unique(total_distance.km)),
                     total_distance.var = var(total_distance.km),
                     max_distance.km = max(max_distance.km),
                     max_distance.var = var(max_distance.km)) %>%
    ungroup() %>%
    filter(duration.hrs > 0) %>%
    filter(!(transit_time.pct > 100 | rest_time.pct > 100)) %>%
    select(-tripID)
  
  all_trips.list[[i]] <- my_trips
   
}

close(pb)
trips_summary <- do.call("rbind", all_trips.list)
trips_summary %<>% filter(duration.days <= 30 & rest_time.pct > 0 & rest_time.pct < 100)

# Change directory back
setwd("C:/Users/gilli/OneDrive - The University of Liverpool/Liverpool postdoc/ALB_foraging/ALB_rpt")
#setwd("C:/Users/ngillies/OneDrive - The University of Liverpool/Liverpool postdoc/ALB_foraging/ALB_rpt")

save(trips_summary, file = paste0("Data_inputs/", my_species, "_", colony_exp, "_individual_trips_summary.RData"))





## APPENDIX --------------------------------------------------------------------

## Fit an HMM to the data

### Prepare the data for HMM

# Remove any trips < 12 hours
gps_labelled.df_hmm <- gps_labelled.df %>% filter(duration.mins >= 720)
gps_labelled.df_hmm %<>% select(longitude, latitude, datetime, ID) %>%
  mutate(ID = as.character(ID)) %>% data.frame()

# Process data for HMM
hmm_data <- prepData(gps_labelled.df_hmm,
                     type = "LL", # longs and lats
                     coordNames = c("longitude", "latitude")) 
head(hmm_data)

# Remove step lengths > 25 - unrealistic speed 
hmm_data <- hmm_data %>%
  filter(!step > 25) 

# Remove NA angles
hmm_data <- hmm_data %>%
  filter(!is.na(angle))

### NOTE: Initial values taken from Clay et al. 2020, J. Anim. Ecol. - need to test new ones

# Assign step lengths 
shape_0 <- c(12.46, 3.95, 0.34)
scale_0 <- c(3.734, 4.44, 0.19) 

# Set zero values to small numbers
ind_zero <- which(hmm_data$step == 0)
if (length(ind_zero) > 0) {
  hmm_data$step[ind_zero] <- runif(length(ind_zero)) / 10000
}
ind_zero <- which(hmm_data$step == "NA")
if (length(ind_zero) > 0) {
  hmm_data$step[ind_zero] <- runif(length(ind_zero)) / 10000
}

stepPar0 <- c(shape_0,scale_0)

# Assign turning angles 
location_0 <- c(0.0033,-0.016, 0.03)
concentration_0 <- c(47.15,  1.16,  39.00)

anglePar0 <- c(location_0, concentration_0)


#### Fit the HMM

stateNames <- c("travel", "search", "rest")

mod_hmm <- fitHMM(
  data = hmm_data,
  nbStates = 3,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(step = stepPar0, angle = anglePar0),
  estAngleMean = list(angle = TRUE),
  stateNames = stateNames
)

#plotPR(mod_hmm) 

# Store model as an .rdata object so you don't have to run from scratch each time
file.out <- paste0("Data_outputs/", my_species, "_", colony_exp, "_HMM.RData")
save(mod_hmm, file = file.out)

#load(file.out)



### Assign behavioural states

# Compute most probable states using viterbi algorithm
hmm_data_out <- mod_hmm$data
hmm_data_out$State <- viterbi(mod_hmm)

# Assign behaviours

## Assess step/angle distributions
ggplot(aes(x = step, fill = as.factor(State)), data = hmm_data_out) + geom_histogram(alpha = 0.5)
ggplot(aes(x = angle, fill = as.factor(State)), data = hmm_data_out) + geom_histogram(alpha = 0.5)

## Label each state and check classification
hmm_data_out$State[hmm_data_out$State == 1] <- "Travel"
hmm_data_out$State[hmm_data_out$State == 2] <- "Search"
hmm_data_out$State[hmm_data_out$State == 3] <- "Rest"

## Check states
table(hmm_data_out$State)

# Calculate percentage time spent in each state 
hmm_data_out %>%
  group_by(State) %>%
  summarize(counts = n()) %>%
  mutate(per = counts / sum(counts) * 100) %>%
  collect()

## Output the data

# Make a counter variable (used in processing in script 2)
hmm_data_out$ID <- as.factor(as.character(hmm_data_out$ID))

hmm_data_out %<>% group_by(ID) %>%
  mutate(counter = row_number(ID)) %>%
  data.frame()

# Rename x and y columns
hmm_data_out <- rename(hmm_data_out, x_lon = x, y_lat = y)

# Remove columns used for processing
cols.rm <- c("ID", "step", "angle")
hmm_data_out[,cols.rm] <- NULL

# Check encoding of variables
hmm_data_out$BirdID <- as.factor(as.character(hmm_data_out$BirdID))
hmm_data_out$TripID <- as.factor(as.character(hmm_data_out$TripID))
factor_vars <- c("Sex", "State")
hmm_data_out[factor_vars] <- lapply(hmm_data_out[factor_vars], factor)
hmm_data_out$DateTime <- as.POSIXct(hmm_data_out$DateTime, format = "%Y-%m-%d %H:%M:%S")

# Output the data
hmm_file <- paste0("Data_inputs/", my_species, "_", colony_exp, "_gps_hmm.RData")
save(hmm_data_out, file = hmm_file)



