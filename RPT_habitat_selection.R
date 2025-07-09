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
colony <- "ker" # cro bi ker
my_species <- "bba" # waal bba

colony_exp <- ifelse(colony == "ker", "kerguelen",
                     ifelse(colony == "cro", "crozet",
                            "birdis"))

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_sst_data.RData"))

# ______________________________ ####

# MODELLING ---------------------------------------------------------------

# * Process data & standardise environment -------

pnts_all %<>%
  mutate(sst = SST - 273.15,
         sst_scaled = scale(sst),
         #bathy_scaled = scale(bathy),
         dist_col_scaled = scale(dist_col),
         used = ifelse(presence == "used", 1, 0)) %>%
  rename(ID = ring) %>%
  select(-SST)

# Remove birds with only a couple of points
pnts_all %<>% 
  group_by(ID, season) %>%
  mutate(n_recs = n(),
         flag = ifelse(n_recs < 20, "flag", "fine")) %>%
  filter(flag == "fine") %>%
  select(-flag)

# Remove NAs and set distance to km
pnts_all %<>% filter(!is.na(sst) & !is.na(used) #& !is.na(bathy)
) %>%
  mutate(dist_col.km = dist_col/1000,
         max_range.km = max_range/1000)

# Make dummy ring_season variable
pnts_all %<>%
  mutate(ID_season = paste0(ID, "_", season)) 


### Plot habitat use --------------------------------------------------------

#### ~# SST # ----------

subset_ids <- sample(unique(pnts_all$ID), 30)

pnts_subset <- pnts_all %>%
  filter(ID %in% subset_ids) %>%
  filter(used == 1) %>%
  group_by(ID) %>%
  mutate(mean_sst = mean(sst)) %>%
  ungroup() %>%
  mutate(ID = factor(ID, levels = unique(ID[order(mean_sst)]))) 

ggplot(pnts_subset, 
       aes(x = sst, y = ID)) +
  geom_density_ridges_gradient(aes(fill = after_stat(x)), scale = 2) +
  scale_fill_viridis_c( name = "SST (°C)" ) +
  labs(x = "Sea surface temperature (°C)",
       caption = "Individual density histograms for use of SST; random sample of 30 individuals.") +
  theme_bw() +
  theme(legend.justification = "top",
        axis.text.y = element_blank(),
        plot.caption = element_text(hjust = 0.5))

#### ~# Bathymetry # ----------

subset_ids <- sample(unique(pnts_all$ID), 30)

pnts_subset <- pnts_all %>%
  filter(ID %in% subset_ids) %>%
  filter(used == 1) %>%
  group_by(ID) %>%
  mutate(mean_bathy = mean(bathy)) %>%
  ungroup() %>%
  mutate(ID = factor(ID, levels = unique(ID[order(mean_bathy)]))) 

ggplot(pnts_subset, 
       aes(x = bathy, y = ID)) +
  geom_density_ridges_gradient(aes(fill = after_stat(x)), scale = 2) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Bathymetry (m)" ) +
  labs(x = "Sea floor depth (m below sea level)",
       caption = "Individual density histograms for use of bathymetry; random sample of 30 individuals.") +
  theme_bw() +
  theme(legend.justification = "top",
        axis.text.y = element_blank(),
        plot.caption = element_text(hjust = 0.5))


#### Variation with distance from colony -------------------------------------

# Bathymetry against maximum range

med_bathy <- pnts_all %>%
  group_by(ID, season, boutID) %>%
  summarise(med_bathy = median(bathy, na.rm = T),
            max_range.km = max_range[1]/1000)

ggplot(med_bathy, aes(y = med_bathy, x = max_range.km)) +
  geom_point() +
  geom_smooth(method = "lm", col = "orange") +
  labs(x = "Maximum range (km)", y = "Median bathymetry (m)",
       caption = "Each point represents an individual trip.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5))


# * Model habitat use -------------------------------------------------------

### Compute models ----------------------------------------------------------

#### ~# SST # ----------

system.time( sst_mod.within <- glmer(
  used ~ sst_scaled + phase + dist_col_scaled + (1|boutID) + # (1|season) +
    (1 + sst_scaled|ID_season),
  data = pnts_all,
  family = binomial(link = "logit") )
)

summary(sst_mod.within)

system.time( sst_mod.between <- glmer(
  used ~ sst_scaled + phase + dist_col_scaled + (1|boutID) + (1|season) +
    (1 + sst_scaled|ID),
  data = pnts_all,
  family = binomial(link = "logit") )
)

summary(sst_mod.between)

save(sst_mod.within, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_sst_scaled_within_glmm.RData"))
save(sst_mod.between, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_sst_scaled_between_glmm.RData"))



#### ~# Bathymetry # ----------
system.time( bathy_mod <- glmer(
  used ~ bathy_scaled + dist_col_scaled + (1|boutID) +
    (1 + bathy_scaled|ID) + (1|season) + (1 + bathy_scaled|ID_season),
  data = pnts_all,
  family = binomial(link = "logit"))
)

relgrad <- with(bathy_mod@optinfo$derivs, solve(Hessian, gradient))
max(abs(relgrad))

summary(bathy_mod)
save(bathy_mod, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_bathy_scaled_glmm.RData"))




# * Interrogate results -----------------------------------------------------

habitat_var <- "sst_scaled" # set to habitat variable: sst, bathy_scaled
fitted_mod.within <- loadRData(paste0("Data_outputs/", my_species, "_", colony_exp, "_", habitat_var, "_within_glmm.RData"))
fitted_mod.between <- loadRData(paste0("Data_outputs/", my_species, "_", colony_exp, "_", habitat_var, "_between_glmm.RData"))

summary(fitted_mod.within); summary(fitted_mod.between)

### Model fit metrics: AUC, sensitivity, specificity ----

# AUC around 0.5 = essential random; closer to 1 = better fit
get_auc(fitted_mod.within); get_auc(fitted_mod.between)

# Check for correlations
car::vif(fitted_mod.within); car::vif(fitted_mod.between)

# Check for random scatter
plot(simulateResiduals(fitted_mod.within))
plot(simulateResiduals(fitted_mod.between))

# Check for patterns and that most points are within 95% confidence intervals
plot_binnedplot(fitted_mod.within)
plot_binnedplot(fitted_mod.between)

# Check ratio is close to 1
overdisp_fun(fitted_mod.within); overdisp_fun(fitted_mod.between)

# Check ratio: >1 may indicate zero inflation (more 0s than expected); < 1 zero deflation
# Small P = evidence to suggest zero-inflation/deflation
testZeroInflation(fitted_mod.within); testZeroInflation(fitted_mod.between)



### Extract & plot random effects ----

## Extract random effects
rnd_eff_within.df <- get_rnd_effects(fitted_mod.within, habitat_var)
rnd_eff_between.df <- get_rnd_effects(fitted_mod.between, habitat_var)

## Slopes & intercepts histograms
par(mfrow = c(2,2))
hist(rnd_eff_within.df$intercept, main = paste0("Within-individiual intercepts: ", habitat_var))
hist(rnd_eff_within.df$slope, main = paste0("Within-individiual slopes: ", habitat_var))
hist(rnd_eff_between.df$intercept, main = paste0("Between-individiual intercepts: ", habitat_var))
hist(rnd_eff_between.df$slope, main = paste0("Between-individiual slopes: ", habitat_var))


# Calculate the predicted probabilities for each ID and habitat value

## Get habitat range
habitat_range <- seq(min(pnts_all[[habitat_var]], na.rm = TRUE), 
                     max(pnts_all[[habitat_var]], na.rm = TRUE), 
                     length.out = 100)

## Make plots
pred_within <- plot_ind_predictions(rnd_eff_within.df, habitat_range, habitat_var, fitted_mod.within)
pred_between <- plot_ind_predictions(rnd_eff_between.df, habitat_range, habitat_var, fitted_mod.between)

### X-lab labeller
x_labeller <- data.frame(var_name = c("bathy_scaled", "sst_scaled"), 
                         label_name = c("Scaled Bathymetry (m)", "ScaledSea surface temperature (°C)"))

pred_within <- pred_within + labs(title = "Between-season", 
                                  x =  x_labeller$label_name[x_labeller$var_name == habitat_var])
pred_between <- pred_between + labs(title = "Within-season",
                                    x = x_labeller$label_name[x_labeller$var_name == habitat_var])

## Output plots
grid.arrange(pred_within, pred_between, ncol = 2)


### Calculate repeatability ----

# Within season #
summary(fitted_mod.within)

# Residual variance
deviance = deviance(fitted_mod.within)
df_resid = df.residual(fitted_mod.within)
var_resid = deviance/df_resid

# Random effect variance
var_comp <- as.data.frame(VarCorr(fitted_mod.within))
var_ID_season.int <- var_comp$vcov[var_comp$grp == "ID_season" & var_comp$var1 == "(Intercept)" & 
                                 is.na(var_comp$var2)]
var_ID_season.slope <- var_comp$vcov[var_comp$grp == "ID_season" & var_comp$var1 == habitat_var &
                                         is.na(var_comp$var2)]
var_ID_season.corr <- var_comp$vcov[var_comp$grp == "ID_season" & var_comp$var1 == "(Intercept)" &
                                       !is.na(var_comp$var2)]
var_boutID <- var_comp$vcov[var_comp$grp == "boutID"]

# Total variance (covariance is multiplied by 2)
var_total = var_ID_season.int + var_ID_season.slope + var_ID_season.corr*2 + var_boutID + var_resid

# Repeatabiltiy estimats
rep_within = var_ID_season.slope / var_total
rep_boutID.within = var_boutID / var_total

rep_within * 100
rep_boutID.within * 100



# Between season #
summary(fitted_mod.between)

# Residual variance
deviance = deviance(fitted_mod.between)
df_resid = df.residual(fitted_mod.between)
var_resid = deviance/df_resid

# Random effect variance
var_comp <- as.data.frame(VarCorr(fitted_mod.between))
var_ID.int <- var_comp$vcov[1]
var_ID.slope <- var_comp$vcov[2]
var_ID.corr <- var_comp$vcov[3]
var_boutID <- var_comp$vcov[4]
var_season <- var_comp$vcov[5]

# Total variance (covariance is multiplied by 2)
var_total = var_ID_season.int + var_ID_season.slope + var_ID_season.corr*2 + var_boutID + var_season + var_resid

# Repeatabiltiy estimats
rep_between = var_ID.slope / var_total
rep_boutID.between = var_boutID / var_total
rep_season = var_season / var_total 

rep_between * 100
rep_boutID.between * 100
rep_season*100







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


# Simulation approach (not working) ---------------------------------------

# Create a raster of a random environmental variable
env <- raster(nrows=100, ncols=100, xmn=col_coords$longitude-2, xmx=col_coords$longitude+2, ymn=col_coords$latitude-1,ymx=col_coords$latitude+2)
env[] <- runif(10000, -80, 180)  
env <- focal(env, w=matrix(1, 5, 5), mean)  
env <- focal(env, w=matrix(1, 5, 5), mean)  

# Generate random birds
## Define class for individual birds
ind <- setClass("ind", slots=c(x="numeric", y="numeric", opt="numeric"))  

## Generate birds with differing habitat preference
n_birds = 30
n_iter = 100
sim_list <- vector(mode = "list", length = n_birds)

for (i in 1:n_birds) {
  
  my_bird <- ind(x = col_coords$longitude, y = col_coords$latitude, opt = runif(1, 0, 100)) 
  
  path_sim <- simulate_movement(my_bird, env, n_iter, 0.1, col_coords$latitude+1.5, 0, col_coords$longitude+1.5, 0)  
  path_sim$ring <- sprintf("B%03d", i)
  colnames(path_sim)[1:2] <- c("longitude", "latitude")
  
  sim_list[[i]] <- path_sim
}

sim_birds <- do.call("rbind", sim_list)

env_rast.spdf <- as(env, "SpatialPixelsDataFrame")
env_ras.df <- as.data.frame(env_rast.spdf)

ggplot() +
  geom_tile(data = env_ras.df, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() +
  geom_path(data = sim_birds, aes(x = longitude, y = latitude, col = ring)) +
  theme_minimal() +
  theme(legend.position = "none")


# Randomly sample environments for each bird & model preference
sim_rings <- unique(sim_birds$ring)
n_pnts = 30
model_outputs.list <- vector(mode = "list", length = length(n_pnts))

for (n in 1:n_pnts) {
  
  sampled.list = list()
  
  for (i in 1:length(sim_rings)) {
    
    sim_birds.subset <- sim_birds %>% filter(ring == sim_rings[i])
    pnts_sampled = sampleRandom(env, size = n*nrow(sim_birds.subset))
    pnts_all <- data.frame(ring = sim_rings[i],
                           used = c(rep(1, nrow(sim_birds.subset)),
                                    rep(0, length(pnts_sampled))),
                           env = c(sim_birds.subset$env_value,
                                   pnts_sampled))
    
    sampled.list[[i]] <- pnts_all
    
  }
  
  sample.df <- do.call("rbind", sampled.list)
  sample.df <- na.omit(sample.df)
  
  hab_pref.glmm <- glmer( used ~ env + (1 + env|ring), data = sample.df,
                          family = binomial(link = "logit") )
  
  glmm.summary <- summary(hab_pref.glmm)
  glmm.varcor <- data.frame(VarCorr(hab_pref.glmm))
  
  stats.df <- data.frame(n_pnts = n,
                         beta_env = glmm.summary$coefficients[2,1],
                         p_env = glmm.summary$coefficients[2,4],
                         slope_env = glmm.varcor[2,4])
  
  model_outputs.list[[n]] <- stats.df
}

ratio_comparison.df <- do.call("rbind", model_outputs.list)
save(ratio_comparison.df, file = paste0("Data_outputs/", my_species, "_", colony, "simulated_ratio_comparison.RData"))


# Group into bouts of searching ----
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

