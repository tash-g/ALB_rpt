## ---------------------------
##
## Script name: extract_foraging_trips.R
##
## Purpose of script: Generalised script to extract foraging trips from GLS and GPS data
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-02-02
##
## Email: gilliesne@gmail.com
##
## ---------------------------
##
## Inputs:
##   - {species}_{colony}_gls.RData
##   - {species}_{colony}_gps.RData
##
## Outputs:
##   - {species}_{colony}_gls_labelled.RData
##   - {species}_{colony}_gps_labelled.RData
##
## ---------------------------


# Set-up =======================================================================

# Load custom functions
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "readxl", "sf", "ggplot2", "data.table")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)
select <- dplyr::select


### Set parameters ----

# Must be manually changed
spec_col <- c("bbal_birdis", "bbal_ker", "waal_birdis", "waal_cro")


# ________________ ####
# * LABEL GPS * ================================================================
  
for(i in 1:length(spec_col)) {
  
  # Load GPS
  my_file <- list.files("Data_inputs/Data_OG/", pattern = paste0(spec_col, "_gps"))
  my_gps <- loadRData(paste0("Data_inputs/Data_OG/", my_file)) 
  
  my_gps %<>%
    mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
    mutate(across(c(latitude, longitude), as.numeric)) %>% filter(!is.na(longitude) & !is.na(latitude)) %>%
    arrange(ring,datetime)
  
  
  ### Label GPS data with colony vs trip ----
  
  my_gps$ringYr <- paste0(my_gps$ring, "_", my_gps$season)
  my_rings <- unique(my_gps$ringYr)
  
  gps_labelled.list <- list()
  
  pb <- txtProgressBar(min = 0, max = length(my_rings), style = 3)
  
  for (i in 1:length(my_rings)) {
    # Progress bar
    setTxtProgressBar(pb, i)
    
    mydat <- subset(my_gps, ringYr == my_rings[i]) %>% distinct()
    
    mydat %<>% 
      filter(!is.na(datetime)) %>%
      arrange(datetime) %>%
      select(c(ringYr, ring, season, datetime, longitude, latitude))
    
    ### Catch problem data
    if ( sum(mydat$longitude) == 0 |
         nrow(mydat) < 5 |
         is.na(mydat$ring)[1] ) {
      next }
    
    #### Calculate speed ----
    mytrip.filtered.df <- mydat
    
    mytrip.filtered.df$calc_speed <- calc_speed(mytrip.filtered.df)
    mytrip.filtered.df %<>% filter(calc_speed < 25)
    
    ### Calculate distance from colony ----
    mytrip.filtered.df$dist_col <- apply(mytrip.filtered.df, 1, function(row) {
      point_coords <- c(as.numeric(row["longitude"]), as.numeric(row["latitude"]))
      geosphere::distHaversine(point_coords, col_coords)
    })
    
    
    ### @CHECKSUMS@ ----
    # Ensure trips start and end within colony #
    
    ## Find minimum distance from colony
    min_col_dist <- median(head(sort(mytrip.filtered.df$dist_col), 10))
    
    if(min_col_dist > 10000) next 
    
    ### Interactively filter out bad trips (uncomment to run) ----
    # plot(mytrip.filtered.df$longitude, mytrip.filtered.df$latitude, type = "l",
    #      main = paste0(my_rings[i], "; ", i, " of ", length(my_rings)))
    # points(mytrip.filtered.df$longitude, mytrip.filtered.df$latitude, col = "red", pch = 20)
    # points(col_coords$longitude, col_coords$latitude, col = "green", pc = 15, cex = 2)
    
    # ## Ask for a decision
    # click <- locator(1)  # Wait for a single click
    # 
    # if (!is.null(click)) { 
    #   remove_ids <- c(remove_ids, my_rings[i])
    #   cat("Marked for removal:", my_rings[i], "\n")
    #   next
    # } 
    # 
    # # Output clean data
    # gps_cleaned.list[[i]] <- mytrip.filtered.df
    
    
    ## Determine whether or not at colony based on distance
    mytrip.filtered.df %<>%
      mutate(loc = ifelse(dist_col < 6000, "colony", "trip")) 
    
    mytrip.filtered.df <- mytrip.filtered.df[order(mytrip.filtered.df$datetime),]
    
    mytrip.filtered.df %<>% 
      mutate(bout = rep(1:length(rle(loc)$lengths), rle(loc)$lengths),
             boutID = paste(loc, ceiling(bout/2), sep = "_")) %>%
      group_by(boutID) %>%
      mutate(duration.mins = as.numeric(difftime(max(datetime), min(datetime), units = "mins")),
             duration.hours = duration.mins/60,
             duration.days = duration.hours/24) 
    
    ## Make sure each trip ends within 40km of colony (allows device to fail 30 mins before arrival assuming max speed 80km/h)
    mytrip.filtered.df %<>% group_by(boutID) %>% 
      mutate(ends_home = ifelse(dist_col[n()] > 40000, "N", "Y"))
    
    mytrip.filtered.df %<>% 
      mutate(ends_home = ifelse(loc == "colony", "Y", ends_home))
    
    mytrip.filtered.df %<>% filter(ends_home == "Y") %>% select(-ends_home)
    
    ### Save data to lists to output ----
    gps_labelled.list[[i]] <- mytrip.filtered.df
    
  }
  
  close(pb)
  
  gps_labelled.df <- do.call("rbind", gps_labelled.list)
  save(gps_labelled.df, file = paste0("Data_inputs/", spec_col, "_gps_labelled.RData"))
  
  
  # SECOND FILTER: Manually check for errors in individual trip tracks (uncomment to run) ----
  
  # load(paste0("Data_inputs/", spec_col, "_gps_labelled.RData"))
  # 
  # gps_labelled.df %<>% mutate(tripID = paste0(ringYr, "_", boutID))
  # my_ids <- unique(gps_labelled.df$tripID)
  # my_ids <- my_ids[!grepl("colony", my_ids)]
  # 
  # remove_trips <- c()
  # 
  # for (i in 1:length(my_ids)) {
  # 
  #   my_trip <- subset(gps_labelled.df, tripID == my_ids[i])
  # 
  #   # Interactively filter out bad trips #
  #   plot(my_trip$longitude, my_trip$latitude, type = "l",
  #        main = paste0(my_ids[i], "; ", i, " of ", length(my_ids)))
  #   points(my_trip$longitude, my_trip$latitude, col = "red", pch = 20)
  #   points(col_coords$longitude, col_coords$latitude, col = "green", pc = 15, cex = 2)
  # 
  #   ## Ask for a decision
  #   click <- locator(1)  # Wait for a single click
  # 
  #   if (!is.null(click)) {
  #     remove_trips <- c(remove_trips, my_ids[i])
  #     cat("Marked for removal:", my_ids[i], "\n")
  #     next
  #   }
  # 
  # }
  # 
  # save(remove_trips, file = paste0("Data_outputs/", spec_col, "_trip_errors.RData"))
  # 
  # gps_labelled.df %<>% filter(!tripID %in% remove_trips)
  # save(gps_labelled.df, file = paste0("Data_inputs/", spec_col, "_gps_labelled.RData"))
  
}

# ________________ ####
# * LABEL GLS * ================================================================

# Load GLS data
my_file <- list.files("Data_inputs/Data_OG/", pattern = paste0(spec_col, "_gls"))
my_gls <- loadRData(paste0("Data_inputs/Data_OG/", my_file)) 

my_gls %<>%
  mutate(season = ifelse(as.numeric(format(as.Date(datetime), "%m")) < 9,
                         as.numeric(format(as.Date(datetime), "%Y")),
                         as.numeric(format(as.Date(datetime), "%Y")) + 1 )) %>%
  select(c(ring, season, datetime, immersion))

### Label GLS data with colony vs trip ----

## Get birdyears
my_gls$birdSeason <- paste(my_gls$ring, my_gls$season, sep = "_")
my_birds <- unique(my_gls$birdSeason)
gls_labelled.list <- vector(mode = "list", length = length(my_birds))

pb <- txtProgressBar(min = 0, max = length(my_birds), style = 3)

for (i in 1:length(my_birds)) {
  
  setTxtProgressBar(pb, i)
  
  # Isolate gls data
  my_gls.loop <- subset(my_gls, birdSeason == my_birds[i])
  my_gls.loop %<>% arrange(datetime) %>% distinct()
  
  ## Skip errors
  if (nrow(my_gls.loop) == 0) next
  if (sum(my_gls.loop$immersion) == 0) next
  
  # Label colony vs trip
  ## Look for dry periods > 8 hours
  my_gls.loop %<>% arrange(datetime) %>%
    mutate(dryWet = ifelse(immersion <= 0.0*max(immersion, na.rm = T), "dry", "wet"),
           dur.mins = rep(rle(dryWet)$lengths, rle(dryWet)$lengths)*10,
           loc = ifelse(dryWet == "dry" & dur.mins > 480, "col", "trip"),
           tripID = rep(1:length(rle(loc)$values), rle(loc)$lengths))
  
  
  gls_labelled.list[[i]] <- my_gls.loop
  
}

close(pb)

gls_labelled.df <- do.call("rbind", gls_labelled.list)

save(gls_labelled.df, file = paste0("Data_inputs/", spec_col, "_gls_labelled.RData"))

