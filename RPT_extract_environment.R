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
              "pROC", "ggridges", "rnaturalearth", "DHARMa", "reticulate") 

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

# Set parameters
colony <- "ker" # cro bi ker
my_species <- "bba" # waal bba

colony_exp <- ifelse(colony == "ker", "kerguelen",
                     ifelse(colony == "cro", "crozet",
                            "birdis"))

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_habitat_avail.RData"))

# ______________________________ ####
# ENVIRONMENT DOWNLOAD --------------------------------------------

# Set parameters
seasons = unique(pnts_all$season)
years = unique(year(pnts_all$datetime))
months = unique(month(pnts_all$datetime))
days = unique(day(pnts_all$datetime))
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)


## * SST -----

### Download file from ERA5 -----

# Set UID and API key
cds.key <- "c446649a-11d0-47a2-8265-858bac139054"
wf_set_key(user = "8942e3d3-fd42-4be5-a07d-6894205633ce", key = cds.key)

# Set request parameters
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels",
  product_type   = "reanalysis",
  format = "grib",
  variable = "sea_surface_temperature", #c("10m_u_component_of_wind", "10m_v_component_of_wind", 
  #"sea_surface_temperature", "total_precipitation")
  year = years,
  month = months,
  day = days,
  time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", 
           "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", 
           "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
  # Specify the area (N, W, S, E)
  area = c(N, W, S, E),
  target = paste0("Data_env/", my_species, "_", colony_exp, "_sst.grib") # update when expand to other vars
)

file <- wf_request(user = "8942e3d3-fd42-4be5-a07d-6894205633ce",
                   request = request,
                   transfer = TRUE,
                   path = "Data_environmental",
                   verbose = TRUE)


#### Visualise SST ------------------------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_habitat_avail.RData"))
grib_data <- rast(paste0("Data_env/", my_species, "_", colony_exp, "_sst.grib"))

time_info <- time(grib_data)
years <- format(as.Date(time_info), "%Y")
months <- format(as.Date(time_info), "%m")

# Aggregate data by year 
annual_sst <- tapp(grib_data, index = years, fun = mean, na.rm = TRUE)

## Convert raster to data frame for ggplot
annual_sst.df <- as.data.frame(annual_sst, xy = TRUE)
colnames(annual_sst.df) <- c("Longitude", "Latitude", paste0("Year_", unique(years)))

## Reshape data for ggplot (annual_sst.df format)
annual_sst.long <- annual_sst.df %>%
  pivot_longer(cols = starts_with("Year_"), names_to = "year", values_to = "SST") %>%
  mutate(year = sub("Year_", "", year),
         SST = SST - 273.15)

## Make the plot
ggplot(data = annual_sst.long, aes(x = Longitude, y = Latitude, fill = SST)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c(name = "SST °C") +
  labs(title = "Annual Mean Sea Surface Temperature") +
  theme_bw() +
  theme(legend.justification = "top")


# Aggregate data by month (mean for each month)
monthly_sst <- tapp(grib_data, index = months, fun = mean, na.rm = TRUE)

## Convert raster to data frame for ggplot
monthly_sst.df <- as.data.frame(monthly_sst, xy = TRUE)
colnames(monthly_sst.df) <- c("Longitude", "Latitude", paste0("Month_", unique(months)))

## Reshape data for ggplot (long format)
monthly_sst.long <- monthly_sst.df %>%
  pivot_longer(cols = starts_with("Month_"), names_to = "month", values_to = "SST") %>%
  mutate(month = sub("Month_", "", month),
         SST = SST - 273.15)

## Reorder months
monthly_sst.long$month <- factor(monthly_sst.long$month, levels = c("12", "01"))

## Make the plot
ggplot(data = monthly_sst.long, aes(x = Longitude, y = Latitude, fill = SST)) +
  geom_tile() +
  facet_wrap(~month, labeller = as_labeller(c("01" = "January", "12" = "December"))) +
  scale_fill_viridis_c(name = "SST °C") +
  labs(title = "Monthly Mean Sea Surface Temperature") +
  theme_bw() +
  theme(legend.justification = "top")


#### Extract SST ------------

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

## Combine the points and their timestamps into a single data.frame
pnts_all$rounded_time <- rounded
pnts_all$index <- seq_len(nrow(pnts_all))  # Add a unique index for reassignment later

## Pre-compute the closest layer index for each unique GPS date
GPSdates = unique(rounded) 
layer_indices <- sapply(GPSdates, function(gps_date) which.min(abs(gps_date - time.index)))

## Create a lookup table for layer indices
layer_lookup <- data.frame(rounded_time = GPSdates, layer_index = layer_indices)

## Join layer indices back to points
pnts_with_layers <- merge(pnts_all, layer_lookup, by = "rounded_time")

# Extract all variables in one pass for each layer
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

## Drop the temporary columns
pnts_all$index <- NULL
pnts_all$rounded_time <- NULL

save(pnts_all, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_sst_data.RData"))

# ______________________________ ####
# * BATHMETRY ---------

### Extract bathymetry ------------

# Set spatial extent
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)

bathy <- getNOAA.bathy(lon1 = W, lon2 = E, lat1 = S, lat2 = N, resolution = 4)

## Check how it looks
plot(bathy, image = TRUE)
scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)
points(pnts_all$Longitude, pnts_all$Latitude, col = "red")

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

save(pnts_all, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_sst-bathy_data.RData"))



### Visualise bathymetry ------------

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_tile(data = bathy.df, aes(x = Longitude, y = Latitude, fill = depth)) +
  geom_sf(data = world, fill = "white", color = "white") +
  geom_point(data = pnts_all %>% filter(used == 1), 
             aes(x = Longitude, y = Latitude), colour = "forestgreen", alpha = 0.5) +
  scale_fill_viridis_c(option = "C", name = "Depth (m)") +
  coord_sf(xlim = c(W+0.5, E-1), ylim = c(S+1, N-1)) +
  theme_bw() +
  theme(legend.justification = "top")



# ______________________________ ####
# * CHLOROPHYLL -------------------------------------------------------------

load(paste0("Data_inputs/", my_species, "_", colony_exp, "_sst-bathy_data.RData"))

### Download chlorophyll data from Copernicus ------------

#### Initialise CMT environment  -----
# Copernicus Marine only interfaces with Python, so need to make a wrapper for it here
# https://help.marine.copernicus.eu/en/articles/8638253-how-to-download-data-via-the-copernicus-marine-toolbox-in-r

# Setup environment - initial (only need to do this once)
# virtualenv_create(envname = "CopernicusMarine")
# virtualenv_install("CopernicusMarine", packages = c("copernicusmarine"))
# reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
# 
# cmt <- import("copernicusmarine")
# 
# # Login function to create your configuration file
# cmt$login("ngillies", ".PVJP-a4-LqcTab")

# Re-active the environment
reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
cmt <- import("copernicusmarine")


#### Download data ----

# Need to run this query year by year as there's no way to subset months otherwise

seasons <- unique(pnts_all$season)

for (i in 1:length(seasons)) {
  
  pb <- txtProgressBar(min = 0, max = length(seasons), style = 3)
  setTxtProgressBar(pb, i)
  
  ## Set parameters
  pnts_subset <- subset(pnts_all, season == seasons[i])
  
  S = min(pnts_subset$Latitude, na.rm = T)
  N = max(pnts_subset$Latitude, na.rm = T)
  W = min(pnts_subset$Longitude, na.rm = T)
  E = max(pnts_subset$Longitude, na.rm = T)
  min_time <- min(pnts_subset$datetime)
  max_time <- max(pnts_subset$datetime)
  
  # Adapt the query for R
  cmt$subset( 
    dataset_id="cmems_mod_glo_bgc_my_0.25deg_P1D-m",
    variables=list("chl"),  
    minimum_longitude=W,  
    maximum_longitude=E, 
    minimum_latitude=S,  
    maximum_latitude=N, 
    start_datetime=format(min_time, "%Y-%m-%dT%H:%M:%S"), 
    end_datetime=format(max_time, "%Y-%m-%dT%H:%M:%S"),  
    minimum_depth=0.5057600140571594,  
    maximum_depth=0.5057600140571594,  
    output_directory = "Data_env/",
    output_filename = paste0(my_species, "_", colony, "_chlorophyll_", seasons[i], ".nc")
  )
}

close(pb)

### Extract chlorophyll -----------------------------------------

# Also needs to be done seasonally
pnts_list <- vector(mode = "list", length = length(seasons))

for (i in 1:length(seasons)) {
  
  # Load the file
  nc_data <- rast(paste0("Data_env/", my_species, "_", colony, "_chlorophyll_", seasons[i], ".nc"))
  
  # split into a separate raster set 
  chl_set = which(grepl("chl_depth=0.50576001", names(nc_data)))
  
  # create separate 'SpatRaster' S4 objects for each  component
  chl_vector <- nc_data[[chl_set]]
  
  # extract the times that are available in the raster files
  time.index = time(chl_vector)
  
  # Round the GPS data timestamps to the nearest day
  pnts_subset <- subset(pnts_all, season == seasons[i])
  formatted_datetime <- format(pnts_subset$datetime, "%Y-%m-%d %H:%M:%S")
  datetime <- ymd_hms(formatted_datetime)
  rounded <- round_date(datetime, "day")
  
  ## Combine the points and their timestamps into a single data.frame
  pnts_subset$rounded_time <- rounded
  pnts_subset$index <- seq_len(nrow(pnts_subset))  # Add a unique index for reassignment later
  
  ## Pre-compute the closest layer index for each unique GPS date
  GPSdates = unique(rounded) 
  layer_indices <- sapply(GPSdates, function(gps_date) which.min(abs(gps_date - time.index)))
  
  ## Create a lookup table for layer indices
  layer_lookup <- data.frame(rounded_time = GPSdates, layer_index = layer_indices)
  
  ## Join layer indices back to points
  pnts_with_layers <- merge(pnts_subset, layer_lookup, by = "rounded_time")
  
  # Extract all variables in one pass for each layer
  system.time(
    results <- pbapply::pblapply(unique(pnts_with_layers$layer_index), function(layer_idx) {
      points <- pnts_with_layers[pnts_with_layers$layer_index == layer_idx, ]
      coords <- data.frame(x = points$Longitude, y = points$Latitude)
      
      chlA <- terra::extract(chl_vector[[layer_idx]], coords)[, 2]
      
      data.frame(index = points$index, chlA = chlA)
    })
  )
  
  # Combine results and assign back to pnts_all
  results_combined <- do.call(rbind, results)
  pnts_subset <- merge(pnts_subset, results_combined, by.x = "index", by.y = "index", all.x = TRUE)
  
  pnts_list[[i]] <- pnts_subset

}

pnts_all <- do.call("rbind", pnts_list)

save(pnts_all, file = paste0("Data_inputs/", my_species, "_", colony_exp, "_sst-bathy-chla_data.RData"))


### Visualise chlorophyll ------------------------------------------------------

for (i in 1:length(seasons)) {
  
  # Load the file
  nc_data <- rast(paste0("Data_env/", my_species, "_", colony, "_chlorophyll_", seasons[i], ".nc"))
  
  time_info <- time(nc_data)
  years <- format(as.Date(time_info), "%Y")
  months <- format(as.Date(time_info), "%m")
  
  # Aggregate data by year 
  annual_chlA <- tapp(nc_data, index = years, fun = mean, na.rm = TRUE)
  
  ## Convert raster to data frame for ggplot
  annual_chlA.df <- as.data.frame(annual_chlA, xy = TRUE)
  colnames(annual_chlA.df) <- c("Longitude", "Latitude", paste0("Year_", unique(years)))
  
  ## Reshape data for ggplot (annual_sst.df format)
  annual_chlA.long <- annual_chlA.df %>%
    pivot_longer(cols = starts_with("Year_"), names_to = "year", values_to = "ChlA") %>%
    mutate(year = sub("Year_", "", year))
  
  ## Make the plot
  ggplot(data = annual_chlA.long, aes(x = Longitude, y = Latitude, fill = ChlA)) +
    geom_tile() +
    facet_wrap(~year) +
    scale_fill_viridis_c(name = "Chl-A") +
    labs(title = "Annual Mean Chlorophyll") +
    theme_bw() +
    theme(legend.justification = "top")
  
  # Aggregate data by month (mean for each month)
  monthly_chlA <- app(nc_data, index = months, fun = mean, na.rm = TRUE)
  
  ## Convert raster to data frame for ggplot
  monthly_chlA.df <- as.data.frame(monthly_chlA, xy = TRUE)
  colnames(monthly_chlA.df) <- c("Longitude", "Latitude", paste0("Month_", unique(months)))
  
  ## Reshape data for ggplot (long format)
  monthly_chlA.long <- monthly_chlA.df %>%
    pivot_longer(cols = starts_with("Month_"), names_to = "month", values_to = "ChlA") %>%
    mutate(month = sub("Month_", "", month))
  
  ## Reorder months
  monthly_chlA.long$month <- factor(monthly_chlA.long$month, levels = c("12", "01"))
  
  ## Make the plot
  ggplot(data = monthly_chlA.long, aes(x = Longitude, y = Latitude, fill = ChlA)) +
    geom_tile() +
    facet_wrap(~month, labeller = as_labeller(c("01" = "January", "12" = "December"))) +
    scale_fill_viridis_c(name = "Chl-A") +
    labs(title = "Monthly Mean Chlorophyll") +
    theme_bw() +
    theme(legend.justification = "top")

}