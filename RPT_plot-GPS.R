## ---------------------------
##
## Script name: RPT_plot-gps.R
##
## Purpose of script: Plot the GPS data to look at individual consistency
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-11-11
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("ALB_FOR_functions.R")
source("ALB_data_processing_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "readxl", "sf", "ggplot2", "trip", "data.table",
              "momentuHMM", "rnaturalearthdata", "sp", "rnaturalearth")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], dependencies =T)
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)
select <- dplyr::select

# Set parameters ===============================================================

# Set colony/species
my_colony <- "bi" # cro bi ker
my_species <- "bba" # waal bba

colony_exp <- ifelse(my_colony == "ker", "kerguelen",
                     ifelse(my_colony == "cro", "crozet",
                            "birdis"))


# ______________________________ ####
# ~ DATA PROCESSING ~ ###########################################################

# PROCESS GPS DATA -------------------------------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_gps_labelled_subset.RData"))

## Calculate landings using speed filter ---------------------------------------

# Remove trips with unrealistic speed
gps_labelled.df %<>% filter(calc_speed < 26)

# Deal with duplicates
gps_labelled.df %<>%
  group_by(ID) %>%
  mutate(datetime = lubridate::round_date(datetime, "10 minutes")) %>%
  distinct() %>%
  arrange(datetime)



# GPS PLOTTING ------------------------------------------------------------

# Get colony locations
if(my_colony == "cro" | my_colony == "ker") { colony.df = data.frame(Lon = 51.706972, Lat = -46.358639) }
#if(my_colony == "ker") { colony.df = data.frame(Lon = 70.25, Lat = -49.68333) }
if(my_colony == "bi") { colony.df = data.frame(Lon = -38.05, Lat = -54) }

# Set projections
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj.utm <- paste0("+proj=laea +lat_0=", colony.df$Lat[1], " +lon_0=", colony.df$Lon[1],
                   " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Colony projection (for plotting)
coordinates(colony.df) <- ~Lon+Lat
proj4string(colony.df) <- proj.dec
colony.sp <- spTransform(colony.df, CRS(proj.utm))
colony.df <- as.data.frame(colony.sp)

## Get the base world plot
world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.utm)

## Plot birds with multiple trips in multiple seasons --------------------------

multi_seasons <- gps_labelled.df %>%
  group_by(ring) %>%
  mutate(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1) %>%
  select(-n_seasons)

## Set projection of GPS data
multi_seasons.proj <- multi_seasons
coordinates(multi_seasons.proj) <- ~ longitude+latitude
proj4string(multi_seasons.proj) <- proj.dec

multi_seasons.sp <- spTransform(multi_seasons.proj, CRS(proj.utm))
multi_seasons.df <- as.data.frame(multi_seasons.sp)

# Loop through target birds
target_birds = unique(multi_seasons$ring) 

for (i in 1:length(target_birds)) {
  
  plotDat <- multi_seasons.df %>% filter(ring == target_birds[i]) %>% filter(loc == "trip")
  padding <- (max(plotDat$coords.x1) - min(plotDat$coords.x1)) * 0.2 
  
  xmin <- min(plotDat$coords.x1) - padding
  xmax <- max(plotDat$coords.x1) + padding
  ymin <- min(plotDat$coords.x2) - padding
  ymax <- max(plotDat$coords.x2) + padding
  
  myplot <- ggplot(data = world2) + 
    geom_sf(fill = "cadetblue", colour = "grey") +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
    coord_sf(crs = proj.utm, xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
            label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
    geom_path(aes(x = coords.x1, y = coords.x2, colour = boutID), size = 1, dat = plotDat) +
    scale_colour_viridis_d(name = NULL) +
    labs(title = paste0("Ring: ", target_birds[i])) +
    facet_grid(.~season) + 
    annotate("point", shape = 17, (colony.df$coords.x1 + 100), colony.df$coords.x2, size = 2) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x.top = element_blank(), 
          axis.title.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y.right = element_blank(), 
          axis.title.y.right = element_blank(),
          axis.title.y.left = element_blank(),
          plot.title = element_text(size = 12, face = "italic"),
          legend.position = "none", 
          strip.text = element_text(face = "bold")) 
  
  # Get number of seasons for plot width
  n_seasons = n_distinct(plotDat$season)
  plot_height = 6
  plot_width = 6*n_seasons
  
  png(file = paste0("Figures/GPS_tracks/Ring/", my_species, "_", my_colony, "_", target_birds[i], ".png"),
      width = plot_width, height = plot_height, units = "in", res = 300)
  print(myplot)
  dev.off()

}


## Plot all tracks for each season --------------------------------------------

## Set projection of GPS data
gps_labelled.df.proj <- gps_labelled.df
coordinates(gps_labelled.df.proj) <- ~ longitude+latitude
proj4string(gps_labelled.df.proj) <- proj.dec

gps_labelled.df.sp <- spTransform(gps_labelled.df.proj, CRS(proj.utm))
gps_labelled.df.df <- as.data.frame(gps_labelled.df.sp)

# Loop through each season
target_seasons = unique(gps_labelled.df$season) 

for (i in 1:length(target_seasons)) {
  
  print(target_seasons[i])
  plotDat <- gps_labelled.df.df %>% filter(season == target_seasons[i]) %>% filter(loc == "trip")
  padding <- (max(plotDat$coords.x1) - min(plotDat$coords.x1)) * 0.2 
  
  xmin <- min(plotDat$coords.x1) - padding
  xmax <- max(plotDat$coords.x1) + padding
  ymin <- min(plotDat$coords.x2) - padding
  ymax <- max(plotDat$coords.x2) + padding
  
  myplot <- ggplot(data = world2) + 
    geom_sf(fill = "cadetblue", colour = "grey") +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
    coord_sf(crs = proj.utm, xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
             label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
    geom_path(aes(x = coords.x1, y = coords.x2, colour = ring), size = 1, dat = plotDat) +
    scale_colour_viridis_d(name = NULL) +
    labs(title = paste0("Season: ", target_seasons[i], "; n = ", n_distinct(plotDat$ring))) +
    annotate("point", shape = 17, (colony.df$coords.x1 + 100), colony.df$coords.x2, size = 2) +
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x.top = element_blank(), 
          axis.title.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y.right = element_blank(), 
          axis.title.y.right = element_blank(),
          axis.title.y.left = element_blank(),
          plot.title = element_text(size = 12, face = "italic"),
          legend.position = "none", 
          strip.text = element_text(face = "bold")) 
  
  png(file = paste0("Figures/GPS_tracks/Season/", my_species, "_", my_colony, "_", target_seasons[i], ".png"),
      width = 7, height = 6, units = "in", res = 300)
  print(myplot)
  dev.off()
  
}

