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
#source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2","gridExtra", "tidyr",
              "sf", "ecmwfr", "lubridate", "terra", "raster", "marmap",
              "ggridges", "rnaturalearth", "reticulate") 

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


# ______________________________ ####

# EXTRACT (loop) ---------------------------------------------------------------

# Run the extraction as a big loop

spec_col <- c("bbal_ker", "bbal_birdis", "waal_birdis", "waal_cro")
ker_coords <- data.frame(longitude = 70.2333, latitude = -49.6833)
bi_coords <- data.frame(longitude = -38.05, latitude = -54.00)
cro_coords <- data.frame(longitude = 51.706972, latitude = -46.358639)

for (n in 1:length(spec_col)) {
  
  print(paste0("Processing ", spec_col[n], "..."))
  
  load(paste0("Data_outputs/", spec_col[n], "_available_pnts_population_level.RData"))
  all_pnts %<>% mutate(year = year(date_hourly))
  pnts_all <- all_pnts
  rm(all_pnts); gc()
  
  # Set spatial extent
  S = min(pnts_all$latitude, na.rm = T)
  N = max(pnts_all$latitude, na.rm = T)
  W = min(pnts_all$longitude, na.rm = T)
  E = max(pnts_all$longitude, na.rm = T)
  
  ## * SST -----------------------------------------------------------------------
  
  grib_data <- rast(paste0("D:/Environmental data/", spec_col[n], "_sst.grib"))
  
  # split into a separate raster set for U and V wind components
  SSTs = which(grepl("temperature", names(grib_data)))
  
  # create separate 'SpatRaster' S4 objects for each wind component
  SST.vector <- grib_data[[SSTs]]
  
  # extract the times that are available in the raster files
  time.index = time(SST.vector)
  
  # Round the GPS data timestamps to the nearest hour
  formatted_datetime <- format(pnts_all$date_hourly, "%Y-%m-%d %H:%M:%S")
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
      coords <- data.frame(x = points$longitude, y = points$latitude)
      
      SST <- terra::extract(SST.vector[[layer_idx]], coords)[, 2]
      
      data.frame(index = points$index, SST = SST)
    })
  )
  
  # Combine results and assign back to pnts_all
  results_combined <- do.call(rbind, results)
  pnts_all <- merge(pnts_all, results_combined, by.x = "index", by.y = "index", all.x = TRUE)
  
  ## Drop the temporary columns
  pnts_all %<>% select(-c(index, rounded_time))
  
  ## Make celsius
  pnts_all %<>% mutate(SST = SST - 273.15)
  
  
  ## * BATHYMETRY ----------------------------------------------------------------
  
  options(timeout = 3000)
  bathy <- getNOAA.bathy(lon1 = W, lon2 = E, lat1 = S, lat2 = N, resolution = 1)
  
  ## Check how it looks
  # plot(bathy, image = TRUE)
  # scaleBathy(bathy, deg = 2, x = "bottomleft", inset = 5)
  # points(pnts_all$Longitude, pnts_all$Latitude, col = "red")
  
  # Convert bathy matrix to dataframe
  lat <- as.numeric(colnames(bathy))  
  lon <- as.numeric(rownames(bathy))  
  depth <- as.vector(bathy)       
  
  bathy.df <- expand.grid(longitude = lon, latitude = lat)
  bathy.df$depth <- depth
  
  save(bathy.df, file = paste0("Data_outputs/", spec_col[n], "_bathymetry_data.RData"))
  
  ## Create a raster from bathymetry matrix
  bathy.raster <- rasterFromXYZ(bathy.df, crs = CRS("+proj=longlat +datum=WGS84"))
  
  # Extract bathymetry values for GPS points
  pnts_all$bathy <- raster::extract(bathy.raster, pnts_all[, c("longitude", "latitude")],
                                    method = "bilinear")
  
  # FOR KERGUELEN ONLY # Calculate distance from shelf break
  
  if (grepl("ker", spec_col[n])) {
    
    contour_lines <- rasterToContour(bathy.raster, levels = -1000)
    plot(contour_lines)
    shelf_break_sf <- st_as_sf(contour_lines)

    pnts_all_sf <- st_as_sf(
      pnts_all,
      coords = c("longitude", "latitude"),  
      crs = 4326  )
    pnts_all$dist_shelf_break <- st_distance(pnts_all_sf, shelf_break_sf)
    pnts_all$dist_shelf_break <- as.numeric(pnts_all$dist_shelf_break) / 1000  # in km
    
    # Assign sign
    pnts_all$dist_shelf_break_signed <- ifelse(pnts_all$bathy > -1000,
                                               pnts_all$dist_shelf_break,
                                               -pnts_all$dist_shelf_break)
  }
  
  ## * CHLOROPHYLL ---------------------------------------------------------------
  
  pnts_all$year <- year(pnts_all$date_hourly)
  years <- unique(pnts_all$year)
  years <- years[years > 1992 & years < 2024]
  pnts_list <- vector(mode = "list", length = length(years))
  
  spec_col.short <- case_when(spec_col[n] == "bbal_ker" ~ "bba_ker",
                              spec_col[n] == "bbal_birdis" ~ "bba_bi",
                              spec_col[n] == "waal_birdis" ~ "waal_bi",
                              spec_col[n] == "waal_cro" ~ "waal_cro")
  
  for (i in 1:length(years)) {
    
    # Load the file
    nc_data <- rast(paste0("Data_env/", spec_col.short, "_chlorophyll_", years[i], ".nc"))
    
    # split into a separate raster set 
    chl_set = which(grepl("chl_depth=0.50576001", names(nc_data)))
    
    # create separate 'SpatRaster' S4 objects for each  component
    chl_vector <- nc_data[[chl_set]]
    
    # extract the times that are available in the raster files
    time.index = time(chl_vector)
    
    # Round the GPS data timestamps to the nearest day
    pnts_subset <- subset(pnts_all, year == years[i])
    formatted_datetime <- format(pnts_subset$date_hourly, "%Y-%m-%d %H:%M:%S")
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
        coords <- data.frame(x = points$longitude, y = points$latitude)
        
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
  
  save(pnts_all, file = paste0("Data_inputs/", spec_col[n], "_sst-bathy-chla_data_population_level.RData"))
  
  print(paste0(spec_col[n], " complete! ", length(spec_col)-n, " more to go."))
  
}




# ______________________________ ####



# DOWNLOAD ---------------------------------------------------------------------

spec_col <- "bbal_birdis"  # "bba_birdis" "bbal_ker" "waal_birdis", "waal_cro"

load(paste0("Data_outputs/", spec_col, "_available_pnts_population_level.RData"))

# Set parameters
years = unique(year(pnts_all$date_hourly))
months = unique(month(pnts_all$date_hourly))
days = unique(day(pnts_all$date_hourly))
S = min(pnts_all$Latitude, na.rm = T)
N = max(pnts_all$Latitude, na.rm = T)
W = min(pnts_all$Longitude, na.rm = T)
E = max(pnts_all$Longitude, na.rm = T)

## * SST -----------------------------------------------------------------------

# Set UID and API key
cds.key <- "c446649a-11d0-47a2-8265-858bac139054"
wf_set_key(user = "8942e3d3-fd42-4be5-a07d-6894205633ce", key = cds.key)

# Set request parameters
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels",
  product_type   = "reanalysis",
  format = "grib",
  variable = "sea_surface_temperature",
  year = years,
  month = months,
  day = days,
  time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", 
           "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", 
           "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
  # Specify the area (N, W, S, E)
  area = c(N, W, S, E),
  target = paste0("Data_env/", my_species, "_", colony, "_sst.grib") # update when expand to other vars
)

file <- wf_request(user = "8942e3d3-fd42-4be5-a07d-6894205633ce",
                   request = request,
                   transfer = TRUE,
                   path = "Data_env",
                   verbose = TRUE)


## * CHLOROPHYLL ---------------------------------------------------------------

# Initialise CMT environment
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


# Download data

# Need to run this query year by year as there's no way to subset months otherwise

pnts_all$year <- format(pnts_all$datetime, "%Y")
seasons <- unique(pnts_all$year)
seasons <- seasons[seasons > 1992 & seasons < 2024]

for (i in 1:length(seasons)) {
  
  ## Set parameters
  pnts_subset <- subset(pnts_all, year == seasons[i])
  
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


# ______________________________ ####
# VISUALISE --------------------------------------------------------------------

spec_col <- "waal_cro"  # "bbal_birdis" "bbal_ker" "waal_birdis", "waal_cro"
spec_col2 <- "waal_birdis"  # "bba_birdis" "bba_kerguelen" "waal_birdis", "waal_cro"

load(paste0("Data_inputs/", spec_col, "_sst-bathy-chla_data_population_level.RData"))
pnts_all %<>% filter(used == 1)

# Set parameters
years = unique(pnts_all$year)

S = min(pnts_all$latitude, na.rm = T)
N = max(pnts_all$latitude, na.rm = T)
W = min(pnts_all$longitude, na.rm = T)
E = max(pnts_all$longitude, na.rm = T)

## * SST -----------------------------------------------------------------------

my_extent <- ext(W, E, S, N)
grib_data <- rast(paste0("D:/Environmental data/", spec_col2, "_sst.grib"))
grib_crop <- crop(grib_data, my_extent)

time_info <- time(grib_crop)
years_sst <- format(as.Date(time_info), "%Y")

annual_sst <- tapp(grib_crop, index = years_sst, fun = mean, na.rm = TRUE, cores = 4)

for (i in 1:nlyr(annual_sst)) {
  
  # Get SST data
  r <- annual_sst[[i]] - 273.15
  sst_df <- as.data.frame(r, xy = TRUE)
  colnames(sst_df) <- c("Longitude", "Latitude", "SST")
  
  year_label <- gsub("X", "", names(annual_sst)[i])
  
  # Get GPS data
  gps_df <- pnts_all %>% filter(year == as.numeric(year_label))
  
  # Make plot
  plot_sst <- ggplot() +
    geom_tile(data = sst_df, aes(x = Longitude, y = Latitude, fill = SST)) +
    geom_path(data = gps_df, aes(x = longitude, y = latitude, group = ring), 
              col = "white", alpha = 0.5) +
    scale_fill_viridis_c(name = "SST (°C)") +
    labs(title = paste("Annual SST", year_label)) +
    coord_fixed() +
    theme_bw()
  
  
  png(filename = paste0("Figures/Environment/sst_", spec_col, "_", gsub("X", "", names(annual_sst)[i]), ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(plot_sst)
  dev.off()

}

## Crozet only 

grib_data <- rast("D:/Environmental data/waal_crozet_sst.grib")

time_info <- time(grib_data)
years_sst <- format(as.Date(time_info), "%Y")


for (i in 15:length(unique(years_sst))) {
  
  # Pick a year to process
  year_target <- unique(years_sst)[i]
  year_idx <- which(years_sst == year_target)
  
  # Subset layers for just that year
  sst_year <- grib_data[[year_idx]]
  
  # Crop to extent
  my_extent <- ext(W, E, S, N)  
  sst_year_crop <- crop(sst_year, my_extent)
  
  # Average across time
  sst_annual <- app(sst_year_crop, fun = mean, na.rm = TRUE, cores = 4)
  
  # Get SST data
  sst_df <- as.data.frame(sst_annual - 273.15, xy = TRUE)
  colnames(sst_df) <- c("Longitude", "Latitude", "SST")
  
  year_label <- year_target
  
  # Get GPS data
  gps_df <- pnts_all %>% filter(year == as.numeric(year_target))
  
  # Make plot
  plot_sst <- ggplot() +
    geom_tile(data = sst_df, aes(x = Longitude, y = Latitude, fill = SST)) +
    geom_path(data = gps_df, aes(x = longitude, y = latitude, group = ring), 
              col = "white", alpha = 0.5) +
    scale_fill_viridis_c(name = "SST (°C)") +
    labs(title = paste("Annual SST", year_label)) +
    coord_fixed() +
    theme_bw()
  
  png(filename = paste0("Figures/Environment/sst_waal_cro_", year_target, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(plot_sst)
  dev.off()
  
}


## * BATHYMETRY ----------------------------------------------------------------

load(paste0("Data_outputs/", spec_col2, "_bathymetry_data.RData"))

world <- ne_countries(scale = "medium", returnclass = "sf")

png(filename = paste0("Figures/Environment/bathymetry_", spec_col, ".png"),
    width = 8, height = 6, units = "in", res = 300)
ggplot() +
  geom_tile(data = bathy.df, aes(x = longitude, y = latitude, fill = depth)) +
  geom_sf(data = world, fill = "white", color = "white") +
  geom_path(data = pnts_all, aes(x = longitude, y = latitude, group = ring), 
            colour = "limegreen", alpha = 0.5) +
  scale_fill_viridis_c(option = "C", name = "Depth (m)") +
  coord_sf(xlim = c(W+0.5, E-1), ylim = c(S+1, N-1)) +
  theme_bw() +
  theme(legend.justification = "top")
dev.off()


## * CHLOROPHYLL ---------------------------------------------------------------

chl_spec_col <- "waal_bi" # "bba_bi" "bba_ker" "waal_bi" "waal_cro"

for (i in 1:length(years)) {
  
  # Load the file
  nc_data <- rast(paste0("Data_env/", chl_spec_col, "_chlorophyll_", years[i], ".nc"))
  
  chl_mean <- app(nc_data, fun = mean, na.rm = TRUE) 
  
  chl_df <- as.data.frame(chl_mean, xy = TRUE)
  colnames(chl_df) <- c("Longitude", "Latitude", "ChlA")
  
  # Get GPS data
  gps_df <- pnts_all %>% filter(used == 1 & year == as.numeric(years[i]))
  
  ## Make the plot
  chl_plot <- ggplot() +
    geom_tile(data = chl_df, aes(x = Longitude, y = Latitude, fill = ChlA)) +
    geom_path(data = gps_df, aes(x = longitude, y = latitude, group = ring), 
              colour = "mediumvioletred", alpha = 0.5) +
    scale_fill_gradient(name = "Chl-A", low = "white", high = "darkgreen",
                        na.value = "grey20") +
    labs(title = paste("Annual Mean Chlorophyll", years[i])) +
    coord_fixed() +
    theme_bw()
  
  png(filename = paste0("Figures/Environment/chlorophyll_", spec_col, "_", years[i], ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(chl_plot)
  dev.off()
  
}




# ______________________________ ####
# APPENDIX ---------------------------------------------------------------------

# * SEA SURFACE HEIGHT -------------------------------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_sst-bathy-chla_data_population_level.RData"))

### Download SSH data from Copernicus ------------

#### Initialise CMT environment  -----

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
## NB Need to use interim dataset for years >= 2021

pnts_all$year <- format(pnts_all$datetime, "%Y")
years <- unique(pnts_all$year)
years_1 <- years[years < 2021 & years > 1992]

years_1 <- years_1[!years_1 %in% c(2001, 2003, 2013, 2014)]

for (i in 1:length(years_1)) {
  
  print(paste0("Processing season ", i, " of ", length(years_1), "..."))
  
  ## Set parameters
  pnts_subset <- subset(pnts_all, year == years_1[i])
  
  S = min(pnts_subset$Latitude, na.rm = T)
  N = max(pnts_subset$Latitude, na.rm = T)
  W = min(pnts_subset$Longitude, na.rm = T)
  E = max(pnts_subset$Longitude, na.rm = T)
  min_time <- min(pnts_subset$datetime)
  max_time <- max(pnts_subset$datetime)
  
  # Adapt the query for R
  cmt$subset( 
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    variables=list("zos"),  
    minimum_longitude=W,  
    maximum_longitude=E, 
    minimum_latitude=S,  
    maximum_latitude=N, 
    start_datetime=format(min_time, "%Y-%m-%dT%H:%M:%S"), 
    end_datetime=format(max_time, "%Y-%m-%dT%H:%M:%S"),  
    minimum_depth=0.49402499198913574,
    maximum_depth=0.49402499198913574,
    output_directory = "Data_env/",
    output_filename = paste0(my_species, "_", colony, "_SSH_", years_1[i], ".nc")
  )
}


years_2 <- years[years >= 2021]

for (i in 1:length(years_2)) {
  
  print(paste0("Processing season ", i, " of ", length(years_2), "..."))
  
  ## Set parameters
  pnts_subset <- subset(pnts_all, year == years_2[i])
  
  S = min(pnts_subset$Latitude, na.rm = T)
  N = max(pnts_subset$Latitude, na.rm = T)
  W = min(pnts_subset$Longitude, na.rm = T)
  E = max(pnts_subset$Longitude, na.rm = T)
  min_time <- min(pnts_subset$datetime)
  max_time <- max(pnts_subset$datetime)
  
  # Adapt the query for R
  cmt$subset( 
    dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1D-m",
    variables=list("zos"),  
    minimum_longitude=W,  
    maximum_longitude=E, 
    minimum_latitude=S,  
    maximum_latitude=N, 
    start_datetime=format(min_time, "%Y-%m-%dT%H:%M:%S"), 
    end_datetime=format(max_time, "%Y-%m-%dT%H:%M:%S"),  
    minimum_depth=0.49402499198913574,
    maximum_depth=0.49402499198913574,
    output_directory = "Data_env/",
    output_filename = paste0(my_species, "_", colony, "_SSH_", years_2[i], ".nc")
  )
}


### Extract SSH -----------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_sst-bathy-chla_data_population_level.RData"))

# Also needs to be done seasonally
pnts_list <- list()

for (i in 1:length(years)) {
  
  print(paste0("Processing year ", i, " of ", length(years), "..."))
  
  # Load the file
  nc_data <- rast(paste0("Data_env/", my_species, "_", colony, "_SSH_", years[i], ".nc"))
  
  # SSH data comes in daily, so need to process by day
  pnts_list.daily <- vector(mode = "list", length = length(names(nc_data)))
  
  for (day in 1:length(names(nc_data))) {
    
    my_day <- names(nc_data)[day]
    
    # split into a separate raster set 
    ssh_set = which(grepl(my_day, names(nc_data)))
    ssh_vector <- nc_data[[ssh_set]]
    
    # extract the times that are available in the raster files
    time.index <- time(ssh_vector) # Extract time from original dataset
    
    # Subst GPS data to this day 
    pnts_subset <- subset(pnts_all, year == years[i] & as.Date(date_hourly) == as.Date(time.index))
    pnts_subset$date_hourly <- as.POSIXct(pnts_subset$date_hourly, tz = "UTC")
    
    ## Convert GPS data to a SpatVector
    pnts_vect <- vect(pnts_subset, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    
    ## Extract SSH values at GPS locations
    ssh_values <- extract(ssh_vector, pnts_vect)
    pnts_subset$SSH <- ssh_values[, 2] 
    
    pnts_list[[length(pnts_list) + 1]] <- pnts_subset
    
  }
  
}

pnts_all <- do.call("rbind", pnts_list)

save(pnts_all, file = paste0("Data_inputs/", my_species, "_", colony_exp, "_sst-bathy-chla-ssh_data_population_level.RData"))


### Visualise SSH ------------------------------------------------------

for (i in 1:length(seasons)) {
  
  # Load the file
  nc_data <- rast(paste0("Data_env/", my_species, "_", colony, "_SSH_", seasons[i], ".nc"))
  
  time_info <- time(nc_data)
  years <- format(as.Date(time_info), "%Y")
  months <- format(as.Date(time_info), "%m")
  
  # Aggregate data by year 
  annual_ssh <- tapp(nc_data, index = years, fun = mean, na.rm = TRUE)
  
  ## Convert raster to data frame for ggplot
  annual_ssh.df <- as.data.frame(annual_ssh, xy = TRUE)
  colnames(annual_ssh.df) <- c("Longitude", "Latitude", paste0("Year_", unique(years)))
  
  ## Reshape data for ggplot (annual_sst.df format)
  annual_ssh.long <- annual_ssh.df %>%
    pivot_longer(cols = starts_with("Year_"), names_to = "year", values_to = "SSH") %>%
    mutate(year = sub("Year_", "", year))
  
  ## Make the plot
  ggplot(data = annual_ssh.long, aes(x = Longitude, y = Latitude, fill = SSH)) +
    geom_tile() +
    facet_wrap(~year) +
    scale_fill_viridis_c(name = "SSH") +
    labs(title = "Sea surface height (m)") +
    theme_bw() +
    theme(legend.justification = "top")
  
  # Aggregate data by month (mean for each month)
  monthly_SSH <- app(nc_data, index = months, fun = mean, na.rm = TRUE)
  
  ## Convert raster to data frame for ggplot
  monthly_SSH.df <- as.data.frame(monthly_SSH, xy = TRUE)
  colnames(monthly_SSH.df) <- c("Longitude", "Latitude", paste0("Month_", unique(months)))
  
  ## Reshape data for ggplot (long format)
  monthly_SSH.long <- monthly_SSH.df %>%
    pivot_longer(cols = starts_with("Month_"), names_to = "month", values_to = "SSH") %>%
    mutate(month = sub("Month_", "", month))
  
  ## Reorder months
  monthly_SSH.long$month <- factor(monthly_SSH.long$month, levels = c("12", "01"))
  
  ## Make the plot
  ggplot(data = monthly_SSH.long, aes(x = Longitude, y = Latitude, fill = SSH)) +
    geom_tile() +
    facet_wrap(~month, labeller = as_labeller(c("01" = "January", "12" = "December"))) +
    scale_fill_viridis_c(name = "SSH") +
    labs(title = "Monthly Mean SSH (m)") +
    theme_bw() +
    theme(legend.justification = "top")
  
}