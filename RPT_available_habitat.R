## ---------------------------
##
## Script name: RPT_available_habitate
##
## Purpose of script: Randomly sample available locations
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-12-10
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2", "lme4", "rptR", "gridExtra", "tidyr",
              "momentuHMM", "sf", "ecmwfr", "lubridate", "terra", "raster", "marmap",
              "pROC", "ggridges", "rnaturalearth", "DHARMa") 

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)

# Make sure using dplyr select
select <- dplyr::select



# Load the data -----------------------------------------------------------

# Set parameters
colony <- "bi" # cro bi ker
my_species <- "waal" # waal bba

colony_exp <- ifelse(colony == "ker", "kerguelen",
                     ifelse(colony == "cro", "crozet",
                            "birdis"))

# Load the GPS data
gps_files <- list.files("Data_inputs/", pattern = paste0(my_species, "_", colony_exp, "_gps"))
my_file <- gps_files[grepl("labelled", gps_files)]

my_gps <- loadRData(paste0("Data_inputs/", my_file)) 

my_trip_meta <- my_gps %>% select(c(ring, season, datetime, dist_col, phase, boutID)) %>%
  group_by(ring, season) %>% mutate(max_range = max(dist_col, na.rm = T)) %>% ungroup()

my_gps %<>%
  mutate(ID = paste0(ring, "_", boutID),
         season = ifelse(as.numeric(format(as.Date(datetime), "%m")) < 9,
                         as.numeric(format(as.Date(datetime), "%Y")),
                         as.numeric(format(as.Date(datetime), "%Y")) + 1 ))

# ______________________________ ####
# Isolate foraging behaviour using HMMs-----------------------------------------

### Prepare dataset ------

my_gps %<>% 
  # Check formatting
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  mutate(across(c(latitude, longitude), as.numeric)) %>% 
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  # Filter for trips only
  mutate(ID = paste0(ring, "_", boutID)) %>% 
  filter(loc == "trip") %>%
  group_by(ID) %>%
  mutate(n_recs = n()) %>%
  # Remove short trips
  filter(n_recs > 3) %>% # loses 28 trips 
  select(-n_recs) %>%
  ungroup() 

my_gps.hmm <- my_gps %>%
  select(ID, datetime, longitude, latitude) %>%
  arrange(ID, datetime) %>%
  relocate(ID, longitude, latitude) %>% 
  data.frame() 

gps_hmm <- prepData(my_gps.hmm,
                    type = "LL", 
                    coordNames = c("longitude", "latitude")) 


# Remove erroneous step lengths that can't be larger than 40 m
hist(gps_hmm$step)
nrow(subset(gps_hmm, step > 40))/nrow(gps_hmm)
gps_hmm <- subset(gps_hmm, step < 40)

### Set 0s to very small numbers 
gps_hmm$step[gps_hmm$step == 0 | is.na(gps_hmm$step)] <- runif(sum(gps_hmm$step == 0 | is.na(gps_hmm$step))) / 10000
gps_hmm <- na.omit(gps_hmm)
gps_hmm %<>% filter(!is.na(angle) & !is.na(step))

save(gps_hmm, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_HMMdata.RData"))



### Fit the HMM  ---------

### NOTE: Initial value finding code taken from Clay et al. 2020 - eventually need to optimise for BBAL

## Assign step lengths based on previously identified initial values
shape_0 <- c(12.42, 4.10, 0.33)
scale_0 <- c(3.62, 4.71, 0.17)
stepPar0 <- c(shape_0, scale_0)

## Assigning turning angles based on previously identified initial values
location_0 <- c(0.00302, 0.00343, 0.0291)
concentration_0 <- c(50.79, 1.27, 44.02)
anglePar0 <- c(location_0, concentration_0)

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_HMMdata.RData"))

stateNames <- c("travel","search", "rest")

my_hmm <- fitHMM(
  data = gps_hmm,
  nbStates = 3,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(step = stepPar0, angle = anglePar0),
  estAngleMean = list(angle = TRUE),
  stateNames = stateNames
)

save(my_hmm, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_HMM.RData"))

## Plot pseudo-residuals
# plotPR(my_hmm)


### Assign behaviours ---------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_HMMdata.RData"))
load(paste0("Data_outputs/", my_species, "_", colony_exp, "_HMM.RData"))

hmm_data_out <- my_hmm$data
hmm_data_out$State <- viterbi(my_hmm)

# Assess step/angle distributions
# ggplot(aes(x = step, fill = as.factor(State)), data = hmm_data_out) + geom_histogram(alpha = 0.5)
# ggplot(aes(x = angle, fill = as.factor(State)), data = hmm_data_out) + geom_histogram(alpha = 0.5)

# Label each state and check classification
hmm_data_out$State[hmm_data_out$State == 1] <- "Travel"
hmm_data_out$State[hmm_data_out$State == 2] <- "Search"
hmm_data_out$State[hmm_data_out$State == 3] <- "Rest"

# Check states
# table(hmm_data_out$State)

# Calculate percentage time spent in each state 
# hmm_data_out %>%
#   group_by(State) %>%
#   summarize(counts = n()) %>%
#   mutate(per = counts / sum(counts) * 100) %>%
#   collect()

# Merge with OG data
hmm_states <- hmm_data_out %>% select(c(ID, datetime, step, angle, State))

my_gps %<>% select(ring, season, datetime, longitude, latitude, dist_col, boutID, ID)
my_gps <- merge(my_gps, hmm_states, by = c("ID", "datetime"), all.x =T)
my_gps %<>% filter(!is.na(step) & !is.na(angle))

save(my_gps, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_labelledHMM.RData"))



# ______________________________ ####
# AVAILABLE HABITAT -------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_labelledHMM.RData"))

# Use maximum foraging range within each season

## Calc max range
my_gps %<>%
  group_by(season, ring) %>%
  mutate(max_range = max(dist_col, na.rm = T)) %>%
  ungroup()

## Get used points
used_pnts <- st_as_sf(my_gps %>% filter(State == "Search"), coords = c ("longitude", "latitude"), crs = 4326) 

# * Create random available points for each used point ------
set.seed(123) 

## Convert used_pnts to match with buffer (in metres)
used_pnts <- st_transform(used_pnts, crs = 3395) 

## Run the sampling function
system.time(
  avail_pnts <- used_pnts %>%
    group_by(ring, season) %>%
    mutate(random_points = purrr::map(geometry, ~sample_points(.x, max_range, 2))) %>%
    unnest(random_points)
)

## Back transform the coordinates
used_pnts <- st_transform(used_pnts, crs = 4326)
avail_pnts <- st_transform(avail_pnts, crs = 4326)

avail_pnts$random_points <- st_set_crs(avail_pnts$random_points, 3395)
avail_pnts$random_points <- st_transform(avail_pnts$random_points, crs = 4326)

## Check it worked
coords <- st_coordinates(avail_pnts$random_points)
avail_pnts$Longitude <- coords[, "X"]
avail_pnts$Latitude <- coords[, "Y"]

used_pnts <- st_transform(used_pnts, crs = 4326) 
coords <- st_coordinates(used_pnts$geometry)
used_pnts$Longitude <- coords[, "X"]
used_pnts$Latitude <- coords[, "Y"]

plot(avail_pnts$Latitude, avail_pnts$Longitude)
points(used_pnts$Latitude, used_pnts$Longitude, col = "red")


### Make a final dataframe for further processing ----------
pnts_used <- data.frame(used_pnts[,c("datetime", "Latitude", "Longitude", "ring")])
pnts_used %<>% select(-geometry) %>% mutate(presence = "used")

pnts_avail <- data.frame(avail_pnts[,c("datetime", "Latitude", "Longitude", "ring")])
pnts_avail %<>% select(-geometry) %>% mutate(presence = "avail")

pnts_all <- rbind(pnts_used, pnts_avail)

save(pnts_all, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_habitat_avail.RData"))


# * Ensure points don't overlap land ----------------------------------------

# Get spatial extent of GPS points
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)


### Make a land raster -------
extent <- ext(W, E, S, N)
resolution <- 0.1 
raster_extent <- rast(extent, res = resolution, crs = "EPSG:4326")

## Download land polygons from rnaturalearth
land <- ne_countries(scale = "medium", returnclass = "sf")

## Rasterize the land polygons
land_raster <- rasterize(vect(land), raster_extent, field = 1, background = 0)

## Plot the land mask
#plot(land_raster, col = c("lightblue", "forestgreen"), main = "Land Mask")


#### Filter points based on land raster -----

# Convert to SpatVector
points_vect <- vect(pnts_all, geom = c("Longitude", "Latitude"), crs = crs(raster_extent))

## Extract raster values at the points' locations
pnts_all$on_land <- extract(land_raster, points_vect)[,2] == 1



#### Resample filtered points ------------------------------------------------

pnts_all <- merge(pnts_all, my_trip_meta, by = c("ring", "datetime"))

pnts_on_land <- st_as_sf(pnts_all %>% filter(on_land == "TRUE"), 
                         coords = c ("Longitude", "Latitude"), crs = 4326) 
pnts_on_land <- st_transform(pnts_on_land, crs = 3395) 

## Run the sampling function
system.time(
  pnts_on_land.resample <- pnts_on_land %>%
    group_by(ring, season) %>%
    mutate(random_points = purrr::map(geometry, ~sample_points(.x, max_range, 10))) %>% # increase sampling to ensure gaps filled
    unnest(random_points)
)

## Back-transform the coordinates
pnts_on_land.resample <- st_transform(pnts_on_land.resample, crs = 4326)
pnts_on_land.resample$random_points <- st_set_crs(pnts_on_land.resample$random_points, 3395)
pnts_on_land.resample$random_points <- st_transform(pnts_on_land.resample$random_points, crs = 4326)

## Convert to dataframe
coords <- st_coordinates(pnts_on_land.resample$random_points)
pnts_on_land.resample$Longitude <- coords[, "X"]
pnts_on_land.resample$Latitude <- coords[, "Y"]

resample.df <- data.frame(pnts_on_land.resample[,c("datetime", "Latitude", "Longitude", "ring")])
resample.df %<>% select(-geometry) %>% mutate(presence = "avail")

## Remove any remaining land points

## Convert to SpatVector
resample.df_vect <- vect(resample.df, geom = c("Longitude", "Latitude"), crs = crs(raster_extent))

## Extract raster values at the points' locations
resample.df$on_land <- extract(land_raster, resample.df_vect)[,2] == 1
resample.df %<>% filter(on_land == FALSE)

##### Bind back into the original DF -----

pnts_all %<>% select(-c(season, dist_col, phase, boutID, max_range)) %>%
  filter(on_land == FALSE)

pnts_all <- rbind(pnts_all, resample.df) %>%
  arrange(ring, datetime)

## Ensure data are balanced
pnts_all %<>%
  group_by(ring, datetime) %>%
  mutate(nrows = n()) %>%
  filter(nrows >= 3) %>%
  arrange(ring, datetime, desc(presence)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  select(-c(nrows, on_land))

## Put metadata back in
pnts_all <- merge(pnts_all, my_trip_meta, by = c("ring", "datetime"), all.x = T)

save(pnts_all, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_habitat_avail.RData"))

# * Validate used vs unused point selection (see Trevail appendix) ----------

