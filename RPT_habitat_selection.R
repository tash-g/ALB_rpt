## ---------------------------
##
## Script name: RPT_habitat_selection
##
## Purpose of script: Create habitat selection models in albatrosses
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
source("ALB_FOR_functions.R")
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2", "lme4", "rptR", "gridExtra", "tidyr",
              "momentuHMM", "sf", "ecmwfr", "lubridate", "terra", "raster", "marmap")

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

## Set parameters
colony <- "bi" # cro bi ker
my_species <- "waal" # waal bba

colony_exp <- ifelse(colony == "ker", "kerguelen",
                     ifelse(colony == "cro", "crozet",
                            "birdis"))


## Load the GPS data
gps_files <- list.files("Data_RPT/", pattern = paste0(my_species, "_", colony_exp, "_gps"))
my_file <- gps_files[grepl("labelled", gps_files)]

my_gps <- loadRData(paste0("Data_RPT/", my_file)) 

my_trip_meta <- my_gps %>% select(c(ring, datetime, distances, phase, boutID))

my_gps %<>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  mutate(across(c(latitude, longitude), as.numeric)) %>% filter(!is.na(longitude) & !is.na(latitude)) %>%
  arrange(ring,datetime) %>%
  ungroup() %>%
  select(c(ring, season, datetime, longitude, latitude, dist_next, distances))

# bbal_ker <- read.csv("Data_inputs/bbal_cr_tracks_15interp_env.csv")
# bbal_ker %<>% select(-c(X, rowID, tel_id, track_tel_id, use_avail)) %>%
#   mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
#          BirdID = as.character(BirdID))

# ______________________________ ####
# Isolate foraging behaviour using HMMs-----------------------------------------

### Prepare dataset ------

my_gps %<>% rename(ID = ring) %>% relocate(ID, longitude, latitude)  %>% data.frame()

gps_hmm <- prepData(my_gps,
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

save(gps_hmm, file = paste0("Data_rpt/", my_species, "_", colony_exp, "_HMMdata.RData"))

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

load(paste0("Data_rpt/", my_species, "_", colony_exp, "_HMMdata.RData"))

stateNames <- c("travel","search", "rest")

my_hmm <- fitHMM(
  data = gps_hmm,
  nbStates = 3,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(step = stepPar0, angle = anglePar0),
  estAngleMean = list(angle = TRUE),
  stateNames = stateNames
)

save(my_hmm, file = paste0("Data_RPT/", my_species, "_", colony_exp, "_HMM.RData"))

# Plot pseudo-residuals
plotPR(hmm_bbal)


### Assign behaviours ---------

load(paste0("Data_rpt/", my_species, "_", colony_exp, "_HMMdata.RData"))
load(paste0("Data_RPT/", my_species, "_", colony_exp, "_HMM.RData"))

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
hmm_states <- hmm_data_out %>% select(c(ID, datetime, step, angle, State)) %>%
  rename(ring = ID)

my_gps <- merge(my_gps, hmm_states, by = c("ring", "datetime"), all.x =T)
my_gps %<>% filter(!is.na(step) & !is.na(angle))

save(my_gps, file = paste0("Data_RPT/", my_species, "_", colony_exp, "_labelledHMM.RData"))


# ______________________________ ####
# Compute 'used' and 'unused' habitat -------------------------------------

load(paste0("Data_RPT/", my_species, "_", colony_exp, "_labelledHMM.RData"))

# Use maximum foraging range within each season

## Calc max range
my_gps %<>%
  group_by(season, ring) %>%
  mutate(max_range = max(distances, na.rm = T))%>%
  ungroup()

## Get used points
used_pnts <- st_as_sf(my_gps %>% filter(State == "Search"), coords = c ("longitude", "latitude"), crs = 4326) 

# Create random available points for each used point
set.seed(123) 

# Vectorized random sampling function
sample_points_vectorized <- function(geometry, max_range, n = 2) {
  buffer <- st_buffer(geometry, dist = max_range[1])
  st_sample(buffer, size = n, type = "random")
}

## Convert used_pnts to match with buffer (in metres)
used_pnts <- st_transform(used_pnts, crs = 3395) 

## Run the sampling function 
system.time(
  avail_pnts <- used_pnts %>%
    group_by(ring, season) %>%
    mutate(random_points = purrr::map(geometry, ~sample_points_vectorized(.x, max_range))) %>%
    unnest(random_points)
)

## Back transform the coordinates
used_pnts <- st_transform(used_pnts, crs = 4326)
avail_pnts <- st_transform(avail_pnts, crs = 4326)

avail_pnts$random_points <- st_set_crs(avail_pnts$random_points, 3395)
avail_pnts$random_points <- st_transform(avail_pnts$random_points, crs = 4326)

### Check it worked
coords <- st_coordinates(avail_pnts$random_points)
avail_pnts$Longitude <- coords[, "X"]
avail_pnts$Latitude <- coords[, "Y"]

used_pnts <- st_transform(used_pnts, crs = 4326) 
coords <- st_coordinates(used_pnts$geometry)
used_pnts$Longitude <- coords[, "X"]
used_pnts$Latitude <- coords[, "Y"]

plot(avail_pnts$Latitude, avail_pnts$Longitude)
points(used_pnts$Latitude, used_pnts$Longitude, col = "red")


### Make a final dataframe for further processing -----
pnts_used <- data.frame(used_pnts[,c("datetime", "Latitude", "Longitude", "ring")])
pnts_used %<>% select(-geometry) %>% mutate(presence = "used")

pnts_avail <- data.frame(avail_pnts[,c("datetime", "Latitude", "Longitude", "ring")])
pnts_avail %<>% select(-geometry) %>% mutate(presence = "avail")

pnts_all <- rbind(pnts_used, pnts_avail)

save(pnts_all, file = paste0("Data_RPT/", my_species, "_", colony_exp, "_habitat_avail.RData"))

# FOR LATER: STRATIFY -----------------------------------------------------

# # Example: Remove points on land
# land_mask <- raster("path_to_land_mask.tif") # Land = 1, Water = 0
# combined_data$on_land <- extract(land_mask, st_coordinates(combined_data))
# 
# # Filter out points on land
# combined_data <- combined_data %>% filter(on_land == 0)


### Validate used vs unused point selection (see Trevail appendix) ----------



# ______________________________ ####
# Integrate environmental data --------------------------------------------

load(paste0("Data_RPT/", my_species, "_", colony_exp, "_habitat_avail.RData"))

## SST -----

### Download file from ERA5 -----

# Set UID and API key
cds.key <- "c446649a-11d0-47a2-8265-858bac139054"
wf_set_key(user = "8942e3d3-fd42-4be5-a07d-6894205633ce", key = cds.key)

# Set parameters
years = unique(year(pnts_all$datetime))
months = unique(month(pnts_all$datetime))
days = unique(day(pnts_all$datetime))
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)
vars =  "sea_surface_temperature" #c("10m_u_component_of_wind", "10m_v_component_of_wind", 
          #"sea_surface_temperature", "total_precipitation")

# Set request parameters
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels",
  product_type   = "reanalysis",
  format = "grib",
  variable = vars,
  year = years,
  month = months,
  day = days,
  time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", 
           "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", 
           "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
  # Specify the area (N, W, S, E)
  area = c(N, W, S, E),
  target = paste0("Data_RPT/", my_species, "_", colony_exp, "_sst.grib") # update when expand to other vars
)

file <- wf_request(user = "8942e3d3-fd42-4be5-a07d-6894205633ce",
                   request = request,
                   transfer = TRUE,
                   path = "Data_environmental",
                   verbose = TRUE)

#### Extract variables ------------

load(paste0("Data_RPT/", my_species, "_", colony_exp, "_habitat_avail.RData"))
grib_data <- rast(paste0("Data_RPT/", my_species, "_", colony_exp, "_sst.grib"))

# split into a separate raster set for U and V wind components
#Us = which(grepl("u wind", names(grib_data)))
#Vs = which(grepl("v wind", names(grib_data)))
SSTs = which(grepl("temperature", names(grib_data)))
#precip = which(grepl("precip", names(grib_data)))

# create separate 'SpatRaster' S4 objects for each wind component
#U.vector <- grib_data[[Us]]
#V.vector <- grib_data[[Vs]]
SST.vector <- grib_data[[SSTs]]
#precip.vector <- grib_data[[precip]]

# extract the v and u components for each coordinate
#time(U.vector) == time(V.vector) # just show that the times are identical, as we will use the same index.

# extract the times that are available in the raster files
time.index = time(SST.vector)

# Round the GPS data timestamps to the nearest hour
formatted_datetime <- format(pnts_all$datetime, "%Y-%m-%d %H:%M:%S")
datetime <- ymd_hms(formatted_datetime)
rounded <- round_date(datetime, "hour")

### OLLIE CODE

# each timestamp in turn.
GPSdates = unique(rounded) # find the raster dates
# create U and V variables for the wind data
pnts_all$sst = NA

# loop through the dates 
system.time(
for(i in 1:length(GPSdates)){
  print(i)
  
  which.layer = which.min(abs(GPSdates[i] - time.index))
  tmp = pnts_all[rounded == GPSdates[i],]
  
  SSTs = terra::extract(SST.vector[[which.layer]], data.frame(x = tmp$Longitude, y = tmp$Latitude))[,2]
  
  pnts_all[rounded == GPSdates[i],]$sst = SSTs
}
)


## BATHMETRY ---------

# Set spatial extent
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)

bathy <- getNOAA.bathy(lon1 = W, lon2 = E, lat1 = S, lat2 = N, resolution = 4)

## Check how it looks
plot(bathy, image = TRUE)
scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)


# Convert bathy matrix to dataframe
lat <- as.numeric(colnames(bathy))  
lon <- as.numeric(rownames(bathy))  
depth <- as.vector(bathy)       

bathy.df <- expand.grid(Longitude = lon, Latitude = lat)
bathy.df$depth <- depth

## Create a raster from bathymetry matrix
bathy.raster <- rasterFromXYZ(bathy.df, crs = CRS("+proj=longlat +datum=WGS84"))

# Extract bathymetry values for GPS points
pnts_all$bathy <- raster::extract(bathy.raster, pnts_all[, c("Longitude", "Latitude")],
                                  method = "bilinear")

save(pnts_all, file = paste0("Data_RPT/", my_species, "_", colony_exp, "_sst-bathy_data.RData"))


### Process data & standardise environment -------

load(paste0("Data_RPT/", my_species, "_", colony_exp, "_sst-bathy_data.RData"))

#source("C:/Users/ngillies/OneDrive - The University of Liverpool/Liverpool postdoc/GPS processing/wind_functions.R")

# pnts_all %<>%
#   mutate(sst = sst - 273.15,
#          windDir = calc_windDir(U_wind, V_wind),
#          dirSigned = calc_signedWindDir(windDir),
#          windSp = calc_windSp(U_wind, V_wind)) %>%
#   mutate(SST_scale = scale(SST),
#          windDir_scale = scale(windDir),
#          dirSigned_scale = scale(dirSigned),
#          windSp_scale = scale(windSp),
#          precip_scale = scale(precip))

pnts_all %<>%
  mutate(sst = sst - 273.15,
         sst_scaled = scale(sst),
         bathy_scaled = scale(bathy))

# ______________________________ ####
# Model habitat use -------------------------------------------------------

### Prepare data for models ----------------------------------------------------------

# Merge metadata
pnts_all <- merge(pnts_all, my_trip_meta %>% rename(ID = ring),
                  by = c("ID", "season", "datetime"), all.x = T)

# Set presence to 0 1
pnts_all %<>% 
  mutate(used = ifelse(presence == "used", 1, 0))

# Get season
pnts_all %<>%
  mutate(season = ifelse(as.numeric(format(as.Date(datetime), "%m")) < 9,
                                 as.numeric(format(as.Date(datetime), "%Y")),
                                 as.numeric(format(as.Date(datetime), "%Y")) + 1 ),)

# Remove NAs
pnts_all %<>% filter(!is.na(sst) & !is.na(used) & !is.na(bathy))

# Group into bouts of searching
# # Sort the data by BirdID and datetime to ensure the observations are in order
# pnts_all %<>%
#   arrange(BirdID, presence, datetime)
# 
# # Calculate the time difference between consecutive observations for each BirdID
# pnts_all %<>%
#   group_by(BirdID, presence) %>%
#   mutate(time_diff = as.numeric(difftime(datetime, lag(datetime), units = "mins"))) %>%
#   ungroup()
# 
# # Assign a group ID to observations within 60 minutes of each other
# pnts_all %<>%
#   group_by(BirdID, presence) %>%
#   mutate(group_id = cumsum(ifelse(is.na(time_diff) | time_diff > 30, 1, 0))) %>%
#   ungroup()
# 
# # Summarise SST
# pnts_all.sum <- pnts_all %>%
#   group_by(BirdID, presence, group_id) %>%
#   summarise(track_id = track_id[1],
#             stt_date = datetime[1],
#             season = season[1],
#             mean_SST = mean(SST, na.rm = T)) %>%
#   mutate(used = ifelse(presence == "used", 1, 0),
#          mean_SST.scaled = scale(mean_SST),
#          season = as.factor(season)) %>%
#   filter(!is.na(mean_SST))

### Compute models ----------------------------------------------------------

library(lme4) 
library(pROC) 

# Fit initial SST model - ~17 minutes
system.time( sst_mod <- glmer(
  used ~ sst + phase + (1|boutID) +
    (1 + sst|ID) + (1|season) + (1 + sst|ID:season),
  data = pnts_all,
  family = binomial(link = "logit")
 )
)


# Fit initial bathymetry model 
system.time( bathy_mod <- glmer(
  used ~ bathy + phase + (1|boutID) +
    (1 + bathy|ID) + (1|season) + (1 + bathy|ID:season),
  data = pnts_all,
  family = binomial(link = "logit")
)
)

summary(bathy_mod)
save(bathy_mod, file = paste0("Data_RPT/", my_species, "_", colony_exp, "_bathy_glmm.RData"))




# Interrogate results -----------------------------------------------------

fitted_mod <- sst_mod # set to focal model
habitat_var <- "sst" # set to habitat variable
 
### Model fit metrics: AUC, sensitivity, specificity ----

# model_data <- model.frame(fitted_mod)
# response <- model_data$used

predicted_probs <- predict(fitted_mod, type = "response")
roc_obj <- roc(pnts_all$used, predicted_probs)
auc_value <- auc(roc_obj)
cat("AUC:", auc_value, "\n") # closer to 0.5 means basically random guess; want something close to 1

### Extract & plot random effects ----
random_effects <- ranef(fitted_mod, condVar = TRUE)
ind_slopes <- random_effects$ID[habitat_var]  
ind_intercepts <- random_effects$ID["(Intercept)"]        

## Plot random effects
rnd_eff.df <- data.frame(cbind(ind_intercepts, ind_slopes))
colnames(rnd_eff.df) <- c("intercept", "slope")
rnd_eff.df$ID = rownames(random_effects$ID)

habitat_range <- seq(min(pnts_all[[habitat_var]], na.rm = TRUE), 
                     max(pnts_all[[habitat_var]], na.rm = TRUE), 
                     length.out = 100)


# Calculate the predicted probabilities for each ID and habitat value
predictions <- expand.grid(ID = rnd_eff.df$ID, habitat_val = habitat_range)

## Merge in slope/intercept
predictions <- merge(predictions, rnd_eff.df, by = "ID", all.x = T)
predictions$logit <- predictions$intercept + predictions$slope * predictions$habitat_val
predictions$probability <- 1 / (1 + exp(-predictions$logit))  

## get global
global_intercept <- fixef(fitted_mod)["(Intercept)"]
global_slope <- fixef(fitted_mod)[[habitat_var]]

predictions.global <- data.frame(habitat_val = habitat_range)
predictions.global$logit <- global_intercept + global_slope * predictions.global$habitat_val
predictions.global$probability <- 1 / (1 + exp(-predictions.global$logit))  

## Plot
ggplot() +
  geom_line(data = predictions, aes(x = habitat_val, y = probability, colour = ID), alpha = 0.6) +
  geom_line(data = predictions.global, aes(x = habitat_val, y = probability), col = "black",
  linetype = "dashed", linewidth = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Habitat value",
       y = "Probability of use")


### Calculate repeatability ----
summary(fitted_mod)

# Residual variance
deviance = deviance(fitted_mod)
df_resid = df.residual(fitted_mod)
var_resid = deviance/df_resid

# Random effect variance
var_comp <- as.data.frame(VarCorr(fitted_mod))
var_ID_season <- var_comp$vcov[var_comp$grp == "ID:season" & var_comp$var1 == "(Intercept)" & 
                                is.na(var_comp$var2)]
var_ID_season_SST <- var_comp$vcov[var_comp$grp == "ID:season" & var_comp$var1 == "sst" &
                                     is.na(var_comp$var2)]

var_ID <- var_comp$vcov[var_comp$grp == "ID" & var_comp$var1 == "(Intercept)" &
                          is.na(var_comp$var2)]
var_ID_sst <- var_comp$vcov[var_comp$grp == "ID" & var_comp$var1 == "sst" &
                             is.na(var_comp$var2)]
var_boutID <- var_comp$vcov[var_comp$grp == "boutID"]

# Total variance
var_total = var_ID_season + var_ID_season_SST + var_ID + var_ID_sst + var_boutID + var_resid

# Repeatabiltiy estimats
rep_within = var_ID_season_SST / var_total
rep_between = var_ID_sst / var_total



# ______________________________ ####
# Troubleshooting ---------------------------------------------------------

# Plot some randomised tracks

temp <- subset(pnts_all, BirdID == "200401")
temp_real <- subset(temp, presence == "used")
temp_sim <- subset(temp, presence == "avail")

plot(temp_real$Latitude, temp_real$Longitude)
points(temp_sim$Latitude, temp_sim$Longitude, col = "red")

# Plot temperature
library(sf)
temp_sf <- st_as_sf(pnts_all %>% filter(presence == "used"), coords = c("Longitude", "Latitude"), crs = "EPSG:4326")
temp_sf <- st_transform(temp_sf, crs(grib_data))

temp_sf.avail <- st_as_sf(pnts_all %>% filter(presence == "avail"), coords = c("Longitude", "Latitude"), crs = "EPSG:4326")
temp_sf.avail <- st_transform(temp_sf.avail, crs(grib_data))


plot(grib_data[[2]])
points(temp_sf$geometry)
points(temp_sf.avail$geometry, col = "red")


# ______________________________ ####
# APPENDIX/Code graveyard -------------------------------------------------


# Alternative code for matching habitat data ------------------------------

### ALTERNATIVE

# Combine the points and their timestamps into a single data.frame
pnts_all$rounded_time <- rounded
pnts_all$index <- seq_len(nrow(pnts_all))  # Add a unique index for reassignment later

# Pre-compute the closest layer index for each unique GPS date
GPSdates = unique(rounded) 
layer_indices <- sapply(GPSdates, function(gps_date) which.min(abs(gps_date - time.index)))

# Create a lookup table for layer indices
layer_lookup <- data.frame(rounded_time = GPSdates, layer_index = layer_indices)

# Join layer indices back to points
pnts_with_layers <- merge(pnts_all, layer_lookup, by = "rounded_time")

# Extract all variables in one pass for each layer
# user  system elapsed 
# 37.99   27.77   66.35 
system.time(
  results <- pbapply::pblapply(unique(pnts_with_layers$layer_index), function(layer_idx) {
    points <- pnts_with_layers[pnts_with_layers$layer_index == layer_idx, ]
    coords <- data.frame(x = points$Longitude, y = points$Latitude)
    
    SST <- terra::extract(SST.vector[[layer_idx]], coords)[, 2]
    
    data.frame(index = points$index, SST = SST)
  })
)

# Combine results and assign back to pnts_all
results_combined <- do.call(rbind, results)
pnts_all <- merge(pnts_all, results_combined, by.x = "index", by.y = "index", all.x = TRUE)

# Drop the temporary columns
pnts_all$index <- NULL
pnts_all$rounded_time <- NULL

### Download SST directly from ERDDAP (NOT WORKING) -----

library(rerddap)

# Download data based on start and end dates(must be within 9 years)
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)

min_time = format(min(pnts_all$datetime), "%Y-%m-%dT00:00:00Z")
max_time = format(max(pnts_all$datetime), "%Y-%m-%dT00:00:00Z")

#rerddap::info(datasetid = "drifter_hourly_qc", url = "https://erddap.aoml.noaa.gov/gdp/erddap/")

system.time(
  OISST_dat <- tabledap(
    x = "drifter_hourly_qc",
    url = "https://erddap.aoml.noaa.gov/gdp/erddap/",
    fields = c("longitude", "latitude", "time", "sst"),
    paste0('time>=', min_time),
    paste0('time<=', max_time),
    paste0('latitude>=', S),
    paste0('latitude<=', N),
    paste0('longitude>=', W),
    paste0('longitude<=', E)
  )
)

OISST_dat <- OISST_dat %>% 
  data.frame %>%
  mutate(sst = sst - 273.15)

save(OISST_dat, file = "ERDDAP_SST_dat.RData")

#### Match environmental data to GPS ----------

# Add an index to OISST data for matching
OISST_dat %<>%
  mutate(time_index = row_number())

# Match closest time - VERY SLOW
GPS_data <- pnts_all %>%
  rowwise() %>%
  mutate(
    closest_time_index = which.min(abs(difftime(OISST_dat$time, datetime, units = "secs"))),
    matched_time = OISST_dat$time[closest_time_index]
  ) %>%
  ungroup()


# Function to calculate Euclidean distance
euclidean_distance <- function(lon1, lat1, lon2, lat2) {
  sqrt((lon1 - lon2)^2 + (lat1 - lat2)^2)
}

# Match closest spatial points
GPS_data <- GPS_data %>%
  rowwise() %>%
  mutate(
    closest_point_index = which.min(
      euclidean_distance(OISST_dat$longitude, OISST_dat$latitude, Longitude, Latitude)
    ),
    matched_longitude = OISST_dat$longitude[closest_point_index],
    matched_latitude = OISST_dat$latitude[closest_point_index]
  ) %>%
  ungroup()

