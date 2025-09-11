## ---------------------------
##
## Script name: extract_trip_behaviours
##
## Purpose: Extract behavioural parameters for foraging effort analyses.
##
## ---------------------------
##
## Dependencies:
##   - extract_foraging_trips -> Splits GLS & GPS datasets into foraging trips
##
## Inputs:
##   - {species}_{colony}_gls_labelled.RData
##   - {species}_{colony}_gps_labelled.RData
##   - {BAS/Chize}_demo_complete.RData
##
## Outputs:
##   - {species}_{colony}_individual_trips_summary.RData
##
## ---------------------------


# Set-up =======================================================================

# Functions for GPS processing and plotting
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "readxl", "sf", "ggplot2", "data.table", "tidyr")

# Install packages not yet installed
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], dependencies =T)
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)
select <- dplyr::select


# ________________ ####
# Process GLS data =============================================================

### Cut GLS to breeding season ----

spec_col <- c("bbal_birdis", "bbal_ker", "waal_birdis", "waal_cro")

for (i in 1:length(spec_col)) {
  
  load(paste0("Data_inputs/", spec_col[i], "_gls_labelled.RData"))
  
  gls_labelled.df %<>%
    mutate(ring = toupper(ring),
           ring = gsub(" ", "", ring),
           ringYr = paste(ring, season, sep = "_")) %>%
    filter(!is.na(ring))
  
  
  #### Load meta data ----
  
  meta_source <- ifelse(grepl("birdis", spec_col[i]), "BAS", "Chize")
  all_demo <- loadRData(paste0("Data_inputs/Data_OG/", meta_source, "_demo_complete.RData"))
  
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
  
  ## Manual breeding dates ----
  
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
  
  ## Match in metadata and filter to breeding ----
  
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
  save(gls_labelled.df, file = paste0("Data_inputs/", spec_col[i], "_gls_labelled_subset.RData"))
  
}



# ________________ ####
# Process GPS data & merge with GPS ============================================

spec_col <- c("bbal_birdis", "bbal_ker", "waal_birdis", "waal_cro")

for (i in 1:length(spec_col)) {

  print(paste0("Processing ", spec_col[i], "..."))
    
  load(paste0("Data_inputs/", spec_col[i], "_gps_labelled.RData")) # dataset filtered manually for bad tracks
  load(paste0("Data_outputs/", spec_col[i], "_labelledHMM.RData")) # dataset labelled with behaviours
  
  # Merge the datasets
  my_gps %<>% select(-c(dist_col, boutID, step, angle))
  
  gps_labelled.df <- merge(gps_labelled.df, my_gps, 
                            by = c("ID", "datetime", "ring", "season", "longitude", "latitude"))
  
  ### Calculate landings based on transitions to flight states ----
  
  ## Deal with potential duplicates
  gps_labelled.df %<>%
    group_by(ID) %>%
    mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>%
    distinct() %>%
    arrange(datetime)
  
  ## Relabel landings & transit to match with GLS
  gps_labelled.df %<>% 
    group_by(ID) %>%
    mutate(behav_state = recode(State,
                                "Travel" = "transit",
                                "Rest" = "rest",
                                "Search" = "aloft"), 
           prev_state = lag(behav_state),
           landing = ifelse( ( behav_state == "aloft" | behav_state == "transit" ) & prev_state == "rest",
                            "take-off",
                            ifelse(behav_state == "rest" & ( prev_state == "aloft" | prev_state == "transit" ),
                                   "landing", "none")),
           landing = ifelse(is.na(landing), "none", landing))
  
  
  ## Calculate distance to next point
  gps_labelled.df %<>%
    group_by(ring, boutID) %>%
    mutate(dist_next = gcd.hf(longitude, latitude, lead(longitude), lead(latitude))) %>%
    ungroup()
  
  #### Create full dataset to integrate with GLS ----
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
  
  ### Load GLS data and merge ----
  load(paste0("Data_inputs/", spec_col[i], "_gls_labelled_subset.RData"))
  
  all_birds <- unique(c(unique(gls_labelled.df$ringYr),
                        unique(gps_full$ringYr)))
  all_trips.list <- vector(mode = "list", length = length(all_birds))
  
  species <- strsplit(spec_col[i], "_")[[1]][1]
  colony <- strsplit(spec_col[i], "_")[[1]][2]
  
  pb <- txtProgressBar(min = 0, max = length(all_birds), style = 3)
  
  for (j in 1:length(all_birds)) {
    
    setTxtProgressBar(pb, j)
    
    mygps <- subset(gps_full, ringYr == all_birds[j]) %>%
      mutate(immersion = NA, dur.mins = NA, 
             dryWet = ifelse(behav_state == "rest", "wet", "dry"),
             source = "GPS") %>%
      ungroup() %>%
      select(c(ring, season, longitude, latitude, datetime, ringYr, phase, loc, rs, 
               lay_date, hatch_date, hatch_code, fail_date, 
               n_landings, behav_state, total_distance.km,
               max_distance.km, immersion, dryWet, source)) 
    
    mygls <- subset(gls_labelled.df, ringYr == all_birds[j]) %>%
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
    
    ### Extract trips ----
    
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
    
    ## Summarise behaviours ----
    my_trips %<>%
      group_by(tripID) %>%
      filter(loc == "trip") %>%
      mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
             prev_state = lag(dryWet),
             bout_time = as.numeric(difftime(datetime, lag(datetime), units = "mins"))) %>%
      dplyr::summarise(species = species,
                       colony = colony,
                       ring = ring[1],
                       season = season[1],
                       phase = phase[1],
                       start_date = datetime[1],
                       end_date = max(datetime, na.rm = T),
                       loc = loc[1],
                       # Duration
                       duration.hrs = as.numeric(difftime(end_date, start_date, units = "hours") + 0.17),
                       duration.days = duration.hrs / 24,
                       # Transitions
                       transitions = sum(dryWet != prev_state, na.rm = TRUE),  # Count state transitions
                       landing = sum(dryWet == "dry" & prev_state == "wet", na.rm = TRUE),
                       takeoff = sum(dryWet == "wet" & prev_state == "dry", na.rm = TRUE),
                       # Behaviours
                       transit_time.hrs = sum(bout_time[behav_state == "transit"], na.rm = T) / 60,
                       rest_time.hrs = sum(bout_time[behav_state == "rest"], na.rm = T) / 60,
                       transit_time.pct = transit_time.hrs/duration.hrs * 100,
                       rest_time.pct = rest_time.hrs/duration.hrs * 100,
                       transit_to_rest.ratio = transit_time.hrs/rest_time.hrs) %>%
      ungroup() %>%
      filter(duration.hrs > 0) %>%
      filter(!(transit_time.pct > 100 | rest_time.pct > 100)) %>%
      select(-tripID)
    
    all_trips.list[[j]] <- my_trips
    
  }
  
  close(pb)
  trips_summary <- do.call("rbind", all_trips.list)
  trips_summary %<>% filter(duration.days <= 30 & rest_time.pct > 0 & rest_time.pct < 100)
  
  # Output data
  save(trips_summary, file = paste0("Data_inputs/", spec_col[i], "_individual_trips_summary.RData"))
  
  
}


# ________________ ####
# Get sample sizes =============================================================

### Overall sample sizes ----

spec_col <- c("bbal_birdis", "bbal_ker", "waal_birdis", "waal_cro")
all_birds.list <- list()

for (i in 1:length(spec_col)) {
  
  species <- strsplit(spec_col[i], "_")[[1]][1]
  colony <- strsplit(spec_col[i], "_")[[1]][2]
  
  load(paste0("Data_inputs/", spec_col[i], "_gls_labelled_subset.RData"))
  gls_summary <- gls_labelled.df %>% select(ring, season, phase) %>% distinct() %>% mutate(logger = "GLS")
  
  load(paste0("Data_inputs/", spec_col[i], "_gps_labelled.RData"))
  gps_summary <- gps_labelled.df %>% select(ring, season, phase) %>% distinct() %>% mutate(logger = "GPS")
  
  all_birds <- rbind(gls_summary, gps_summary) %>% mutate(species = species, colony = colony)
  all_birds.list[[i]] <- all_birds
  
}

all_birds.df <- do.call("rbind", all_birds.list)

all_birds.df %>% group_by(species, colony) %>% summarise(n_sample = n_distinct(ring))

by_logger <- all_birds.df %>% select(ring, season, logger) %>% distinct() %>% group_by(ring, season) %>% mutate(logger = ifelse(n_distinct(logger) > 1, "both", logger))
table(by_logger$logger)

by_individual <- by_logger %>% ungroup() %>% select(-season) %>% distinct()
table(by_individual$logger)

### GLS sample sizes ----

spec_col <- c("bbal_birdis", "bbal_ker", "waal_birdis", "waal_cro")
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



### GPS sample sizes ----

sample_gps.list <- list()
runs_gps.list <- list()

for (i in 1:length(spec_col)) {
  
  load(paste0("Data_inputs/", spec_col[i], "_gps_labelled.RData"))
  
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
