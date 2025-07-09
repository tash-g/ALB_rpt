## ---------------------------
##
## Script name: RPT_habitat_selection
##
## Purpose of script: Create habitat selection models in albatrosses
##
## NOTE: Now working entirely in the R Markdown file.
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
packages <- c("dplyr", "magrittr", "ggplot2", "lme4", "ggpubr", "tidyr",
              "sf", "ecmwfr", "ggridges", "rnaturalearth", "DHARMa",
              "geosphere", "purrr", "lmerTest", "sp", "dtw", "brms") 


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


# Set parameters -----
colony <- "bi" # cro bi ker
my_species <- "waal" # waal bba

colony_exp <- ifelse(colony == "ker", "kerguelen",
                     ifelse(colony == "cro", "crozet",
                            "birdis"))


# Set colony coordinates for bearing calculation
ker_coords <- data.frame(longitude = 70.2333, latitude = -49.6833)
bi_coords <- data.frame(longitude = -38.05, latitude = -54.00)
cro_coords <- data.frame(longitude = 51.706972, latitude = -46.358639)

if (colony == "ker") {
  col_coords <- ker_coords
} else if (colony == "bi") {
  col_coords <- bi_coords
} else if (colony == "cro") {
  col_coords <- cro_coords }


# Load the data -----------------------------------------------------------

gps_files <- list.files("Data_inputs/", pattern = paste0(my_species, "_", colony_exp, "_gps"))
my_file <- gps_files[grepl("labelled", gps_files)]

my_gps <- loadRData(paste0("Data_inputs/", my_file)) 

# Sample sizes
my_gps %>% filter(loc == "trip") %>% group_by(ring) %>% mutate(n_trips = n_distinct(tripID)) %>% filter(n_trips > 1) %>% 
  ungroup() %>% group_by(season) %>% summarise(n_birds = n_distinct(ring))

my_trip_meta <- my_gps %>% select(c(ring, season, datetime, dist_col, phase, boutID)) %>%
  group_by(ring, season) %>% mutate(max_range = max(dist_col, na.rm = T)) %>% ungroup()

my_gps %<>%
  mutate(ID = paste0(ring, "_", season, "_", boutID),
         season = ifelse(as.numeric(format(as.Date(datetime), "%m")) < 9,
                         as.numeric(format(as.Date(datetime), "%Y")),
                         as.numeric(format(as.Date(datetime), "%Y")) + 1 ))

## Filter out individuals with only one trip
my_gps %<>% 
  filter(loc == "trip") %>%
  group_by(ring) %>%
  mutate(n_trips = n_distinct(tripID)) %>%
  filter(n_trips > 1) %>%
  select(-n_trips)

my_gps %>% group_by(season) %>% summarise(n_birds = n_distinct(ring))


## Make a metadata dataframe
my_meta <- my_gps %>%
  group_by(ID) %>%
  select(c(ring, season, boutID, duration.mins, duration.hours, ID, phase)) %>%
  distinct()


# MANUAL FILTER:
# Error in one of the BI BBAL tracks
x <- subset(my_gps, tripID == "1301143_2002_trip_5")

if(nrow(x) != 0) {
  
  my_gps %<>% filter(!(tripID == "1301143_2002_trip_5" & 
             datetime >= as.POSIXct("2002-03-05 12:54:00", tz = "GMT") & 
               datetime <= as.POSIXct("2002-03-05 14:54:00", tz = "GMT")))
  
}


# ______________________________ ####

# Plot the GPS tracks -----------------------------------------------------

# Set projections (use same for Ker + Cro)
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

if(colony == "bi") {
  proj.utm = "+proj=laea +lat_0=-54 +lon_0=-38.03 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
} else {
  proj.utm = "+proj=laea +lat_0=-46.358639 +lon_0=51.706972 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
}

# Get the base world plot
world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.utm)

# Project GPS data
my_gps.proj <- my_gps %>% arrange(ringYr, datetime)
coordinates(my_gps.proj) <- ~ longitude+latitude
proj4string(my_gps.proj) <- proj.dec

my_gps.proj.sp <- spTransform(my_gps.proj, CRS(proj.utm))
my_gps.proj.df <- as.data.frame(my_gps.proj.sp)

# Project the colony
coordinates(col_coords) <- ~longitude+latitude
proj4string(col_coords) <- proj.dec
colony.sp <- spTransform(col_coords, CRS(proj.utm))
colony.df <- as.data.frame(colony.sp)


# Build the plot

## Set min/max limits
y_min = min(my_gps.proj.df$coords.x2 - 100000)
y_max = max(my_gps.proj.df$coords.x2 + 100000)
x_min = min(my_gps.proj.df$coords.x1 - 100000)
x_max = max(my_gps.proj.df$coords.x1 + 100000)
 
# y_min = y_min/4
# y_max = y_max/4
# x_min = x_min/4
# x_max = x_max/4

gps_tracks.plot <- ggplot(data = world2) + 
  geom_sf(fill = "cadetblue", colour = "cadetblue") +
  coord_sf(crs = proj.utm, xlim = c(x_min, x_max), ylim = c(y_min, y_max),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  geom_path(aes(x = coords.x1 , y = coords.x2, col = ringYr), size = 0.8,
            dat = my_gps.proj.df, group = "identity") +
  scale_x_continuous(breaks = seq(-180, 180, by = 5)) +
  scale_y_continuous(breaks = seq(0, -90, by = -5)) + 
  scale_colour_viridis_d() +
  annotate("point", shape = 17, size = 3, (colony.df$coords.x1 + 100), colony.df$coords.x2) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y.left = element_text(size = 10),
        axis.text.x.top = element_text(size = 10),
        legend.position = "none")


# ______________________________ ####
# Approach 1: Bearings from colony -----------------------------------------------

### Calculate bearing differences from colony to distal point ------

## Nest the data by ID
nested_trips <- my_gps %>%
  filter(loc == "trip") %>%
  group_by(ID) %>%
  nest()

## Apply bearing function to each trip
trip_bearings <- nested_trips %>%
  mutate(bearing_col = map_dbl(data, ~ bearing_from_col(.x, col_coords$longitude, col_coords$latitude))) %>%
  select(ID, bearing_col) 

## Put metadata back in
trip_bearings <- merge(trip_bearings, my_meta, by = "ID", all.x = T)

## Convert bearings to radians
trip_bearings %<>% mutate(bearing_rad = deg2rad(bearing_col))


# Within-individual bearings #

# Make sure boutIDs are unique to season
trip_bearings %<>% mutate(boutID = paste0(season, "_", boutID))

# # Compute differences and remove any that occur in different breeding phases
# bearings_within <- trip_bearings %>%
#   group_by(ring) %>%
#   summarise(trip_comparisons = list(combn(boutID, 2, paste, collapse = " vs ")),  # Store trip comparisons
#             bearing_diff = list(combn(bearing_col, 2, function(x) angle_diff(x[1], x[2]))),
#             .groups = "drop") %>%
#   unnest(cols = c(trip_comparisons, bearing_diff)) %>%
#   separate(trip_comparisons, into = c("tripID_1", "tripID_2"), sep = " vs ") %>%
#   left_join(trip_bearings %>% select(ring, boutID, season, phase), 
#             by = c("ring", "tripID_1" = "boutID")) %>%
#   rename(season_1 = season, phase_1 = phase) %>%
#   left_join(trip_bearings %>% select(ring, boutID, season, phase), 
#             by = c("ring", "tripID_2" = "boutID")) %>%
#   rename(season_2 = season, phase_2 = phase) %>%
#   mutate(comp_ind = "within_ring",
#          comp_season = ifelse(season_1 == season_2, "within_season", "between_season"),
#          comp_phase = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")) %>%
#   arrange(ring, season_1, season_2) %>%
#   filter(comp_phase == "same_phase") %>%
#   select(-c(phase_1, phase_2, comp_phase)) %>%
#   # Add ring 2 for compatability with larger df
#   mutate(ring_2 = ring) %>%
#   rename(ring_1 = ring) %>%
#   relocate(ring_1, tripID_1, season_1, ring_2, tripID_2, season_2, bearing_diff, comp_ind, comp_season) %>%
#   # Remove duplicate comparisons
#   mutate(tripID_1.tmp = paste0(ring_1, tripID_1),
#          tripID_2.tmp = paste0(ring_2, tripID_2),
#          comp_ID = paste0(pmin(tripID_1.tmp, tripID_2.tmp), "_", pmax(tripID_1.tmp, tripID_2.tmp))) %>%
#   distinct(comp_ID, .keep_all = TRUE) %>%
#   select(-c(comp_ID, tripID_1.tmp, tripID_2.tmp)) 
# 
# 
# 
# # Among-individual differences #
# 
# # Make the dataframe of possible comparisons to start
# bearings_between.df <- expand.grid(
#   ID1 = unique(trip_bearings$ID),
#   ID2 = unique(trip_bearings$ID)) %>%
#   separate(ID1, into = c("ring_1", "season_1", "dummy1", "tripID_1"), sep = "_", remove = F) %>%
#   separate(ID2, into = c("ring_2", "season_2", "dummy2", "tripID_2"), sep = "_", remove = F) %>%
#   select(-c(dummy1, dummy2)) %>%
#   mutate(tripID_1 = paste0(season_1, "_trip_", tripID_1),
#          tripID_2 = paste0(season_2, "_trip_", tripID_2),
#          # Remove duplicate comparisons
#          tripID_1.tmp = paste0(ring_1, tripID_1),
#          tripID_2.tmp = paste0(ring_2, tripID_2),
#          comp_ID = paste0(pmin(tripID_1.tmp, tripID_2.tmp), "_", pmax(tripID_1.tmp, tripID_2.tmp))) %>%
#   distinct(comp_ID, .keep_all = TRUE) %>%
#   select(-c(comp_ID, tripID_1.tmp, tripID_2.tmp)) %>%
#   filter(ring_1 != ring_2) %>%
#   # Get phase back and filter for same breeding phase
#   left_join(trip_bearings %>% select(ring, boutID, phase, bearing_col) %>% 
#               rename(ring_1 = ring, tripID_1 = boutID),
#             by = c("ring_1", "tripID_1")) %>%
#   rename(phase_1 = phase, bearing_1 = bearing_col) %>%
#   left_join(trip_bearings %>% select(ring, boutID, phase, bearing_col) %>% 
#               rename(ring_2 = ring, tripID_2 = boutID),
#             by = c("ring_2", "tripID_2")) %>%
#   rename(phase_2 = phase, bearing_2 = bearing_col) %>%
#   mutate(phase_comp = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")) %>%
#   filter(phase_comp == "same_phase") %>%
#   select(-c(phase_comp, phase_1, phase_2))
# 
# 
# ## Downsample the comparison data and calculate the differences
# set.seed(817)
# n = 10
# 
# bearings_between <- bearings_between.df %>% 
#   sample_n(size = n * nrow(bearings_within), replace = FALSE)
# 
# ## Calculate comparisons
# bearings_between %<>% 
#   mutate(bearing_diff = angle_diff(bearing_1, bearing_2),
#          comp_ind = "between_ring",
#          comp_season = ifelse(season_1 == season_2, "within_season", "between_season")) %>%
#   select(c(ring_1, tripID_1, season_1, ring_2, tripID_2, season_2, bearing_diff, comp_ind, comp_season)) %>%
#   relocate(ring_1, tripID_1, season_1, ring_2, tripID_2, season_2, bearing_diff, comp_ind, comp_season)
# 
# # Combine results #
# bearing_differences <- rbind(bearings_within, bearings_between)
# 
# save(bearing_differences, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_bearing_differences.RData"))



### Plot bearings in circular plot -----------------------------------

trip_bearings %<>% mutate(bearing_360 = ifelse(bearing_col < 0, bearing_col + 360, bearing_col))

breaks = seq(0, 360, 60)

bearings.plot <- ggplot(trip_bearings, aes(x = bearing_360)) +
  geom_histogram(binwidth = 10,
                 fill = "skyblue", color = "black", alpha = 0.7) + 
  coord_polar(start = 0) +
  scale_x_continuous("", limits = c(0, 360), 
                     breaks = breaks, #c(0, 90, 180, 270), # breaks,
                     labels = paste0(breaks, "\u00B0")) +#c("N", "E", "S", "W")) + # paste0(breaks, "\u00B0")) +
  labs(y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks.y = element_line(color = "black"))



### Repeatability of bearings ------------------------------------

# Integrate meta
trip_bearings.rpt <- merge(trip_bearings, all_ind.meta %>% filter(species == toupper(my_species)), by = "ring")

# Fit the model
bearings_rpt.brms <- brm(
  formula = bf(bearing_rad ~ (1|ring), family = von_mises()),
  data = trip_bearings.rpt,
  chains = 4, iter = 2000, cores = 4)

summary(bearings_rpt.brms)

save(bearings_rpt.brms, file = paste0("Data_outputs/bearings_rpt_brms_", colony, "_", my_species, ".RData"))

posterior_samples <- posterior_samples(bearings_rpt.brms, pars = c("sd_ring__Intercept", "kappa"))

sd_individual <- posterior_samples$sd_ring__Intercept  
kappa <- posterior_samples$kappa  

# Approximate within-individual variance - kappa is the concentration parameter; higher K = less variance (and 
# vice versa) so we divide it by 1 as it is inversely related to dispersion 
var_within <- 1 / kappa  

# Compute repeatability
repeatability <- sd_individual^2 / (sd_individual^2 + var_within)

# Summary statistics
med_R.bearing <- median(repeatability) 
ci_R.bearing <- quantile(repeatability, c(0.025, 0.975))  



### Model bearing differences --------------------------------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_bearing_differences.RData"))

## Merge sex data
load("Data_inputs/ind_meta.RData")

bearing_differences.sex <- merge(bearing_differences, all_ind.meta %>% select(-c(birthYr, species)) %>% rename(ring_1 = ring), 
                by = "ring_1", all.x = T)


##### Within-individual, across seasons ---------

# Fit the model
bearings_within.lmm <- lmer(bearing_diff ~ comp_season + (1+comp_season|ring_1), 
                            data = bearing_differences %>% filter(comp_ind == "within_ring"))

summary(bearings_within.lmm)


# Summarise means
bearings_within.summary <- summary(bearings_within.lmm)

bearings_within.means <- data.frame(comp_season = c("between_season", "within_season"),
                                    mean = c(bearings_within.summary$coefficients[1,1], bearings_within.summary$coefficients[1,1] + bearings_within.summary$coefficients[2,1]),
                                    sd = c(bearings_within.summary$coefficients[1,2], bearings_within.summary$coefficients[2,2]))



# Plot #
bearings_within.plot <- ggplot() +
  geom_histogram(data = bearing_differences %>% filter(comp_ind == "within_ring"), aes(x = abs(bearing_diff))) +
  geom_segment(data = bearings_within.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = bearings_within.means, aes(x = mean+sd, xend = mean+sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = bearings_within.means, aes(x = mean-sd, xend = mean-sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_season, ncol = 1, labeller = labeller(comp_season = c(
    "within_season" = "Within Season",
    "between_season" = "Between Season" ))) +
  labs(y = "Frequency", x = "Angular bearing differences, \u00B0",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear mixed effects model.") +
  theme_bw()


# Calculate individual repeatability
bearings_within.varcor <- as.data.frame(VarCorr(bearings_within.lmm))

ind_rpt.intercept <- bearings_within.varcor$vcov[1] / (bearings_within.varcor$vcov[1] + bearings_within.varcor$vcov[4])
ind_rpt.slope <- bearings_within.varcor$vcov[2] / (bearings_within.varcor$vcov[2] + bearings_within.varcor$vcov[4])


within_individual.random <- ranef(bearings_within.lmm)$ring_1

within_individual.slopes <- within_individual.random$comp_seasonwithin_season
within_individual.intercept <- within_individual.random$`(Intercept)`


# Make a summary output
bearings_within.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_within.lmm"),
  intercept = bearings_within.summary$coefficients[1,1],
  intercept_sd = bearings_within.summary$coefficients[1,2],
  t(bearings_within.summary$coefficients[2,]),
  intercept_var = bearings_within.varcor$vcov[1],
  slope_var = bearings_within.varcor$vcov[2],
  resid_var = bearings_within.varcor$vcov[4]
) )


bearings_within.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))


##### Between-individual, across seasons  -------------------------------

# Fit the model
bearings_between.lm <- lm(bearing_diff ~ comp_ind, data = bearing_differences)
summary(bearings_between.lm)

## Summarise mean differences
bearings_between.summary <- summary(bearings_between.lm)
bearings_between.means <- data.frame(comp_ind = c("between_ring", "within_ring"),
                                     mean = c(bearings_between.summary$coefficients[1,1], bearings_between.summary$coefficients[1,1] + bearings_between.summary$coefficients[2,1]),
                                     sd = c(bearings_between.summary$coefficients[1,2], bearings_between.summary$coefficients[2,2]))

# Plot #
bearings_between.plot <- ggplot() +
  geom_histogram(data = bearing_differences, aes(x = bearing_diff)) +
  geom_segment(data = bearings_between.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = bearings_between.means, aes(x = mean+sd, xend = mean+sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = bearings_between.means, aes(x = mean-sd, xend = mean-sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_ind, ncol = 1, labeller = labeller(comp_ind = c(
    "within_ring" = "Within individual",
    "between_ring" = "Between individual" ))) +
  labs(y = "Frequency", x = "Angular bearing differences, \u00B0",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear model.")  +
  theme_bw()


# Make summary output
bearings_between.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between.lm"),
  intercept = bearings_between.summary$coefficients[1,1],
  intercept_sd = bearings_between.summary$coefficients[1,2],
  t(bearings_between.summary$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )


bearings_between.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))

##### Between-individual, within seasons  -------------------------------

# Need to ensure that only calculate when have multiple birds per season
bearing_differences.filter <- bearing_differences %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")

# Fit the model
bearings_between_withinseason.lmm <- lmer(bearing_diff ~ comp_ind + (1|season_1), data = bearing_differences.filter)
summary(bearings_between_withinseason.lmm)

## Extract mean differences from lmm
bearings_between_withinseason.summary <- summary(bearings_between_withinseason.lmm)

bearings_between_withinseason.means <- data.frame(comp_ind = c("between_ring", "within_ring"),
                                            mean = c(bearings_between_withinseason.summary$coefficients[1,1], bearings_between_withinseason.summary$coefficients[1,1] + bearings_between_withinseason.summary$coefficients[2,1]),
                                            sd = c(bearings_between_withinseason.summary$coefficients[1,2], bearings_between_withinseason.summary$coefficients[2,2]))

# Plot #
bearings_between_withinseason.plot <- ggplot() +
  geom_histogram(data = bearing_differences %>% filter(comp_season == "within_season"), aes(x = abs(bearing_diff))) +
  geom_segment(data = bearings_between_withinseason.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = bearings_between_withinseason.means, aes(x = mean+sd, xend = mean+sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = bearings_between_withinseason.means, aes(x = mean-sd, xend = mean-sd, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_ind, ncol = 1, labeller = labeller(comp_ind = c(
    "within_ring" = "Within individual",
    "between_ring" = "Between individual" ))) +
  labs(y = "Frequency", x = "Angular bearing differences, \u00B0",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear mixed effects model.") +
  theme_bw()


# Calculate individual repeatability
bearings_between_withinseason.varcor <- as.data.frame(VarCorr(bearings_between_withinseason.lmm))

ind_rpt.intercept <- bearings_between_withinseason.varcor$vcov[1] / (bearings_between_withinseason.varcor$vcov[1] + bearings_between_withinseason.varcor$vcov[2])


# Make summary output
bearings_between_withinseason.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between_withinseason.lmm"),
  intercept = bearings_between_withinseason.summary$coefficients[1,1],
  intercept_sd = bearings_between_withinseason.summary$coefficients[1,2],
  t(bearings_between_withinseason.summary$coefficients[2,]),
  intercept_var = bearings_between_withinseason.varcor$vcov[1],
  slope_var = NA,
  resid_var = bearings_between_withinseason.varcor$vcov[2]
) )


bearings_between_withinseason.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))



### By sex ------------------------------------------------------------------

bearing_differences.F = bearing_differences.sex %>% filter(sex == "F")
bearing_differences.M = bearing_differences.sex %>% filter(sex == "M")

# Need to ensure that only calculate when have multiple birds per season
bearing_differences.filter.F <- bearing_differences.F %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")

bearing_differences.filter.M <- bearing_differences.M %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")


###### Between-individual, across seasons  -------------------------------

## FEMALES ##
bearings_between.lm.F <- lm(bearing_diff ~ comp_ind, data = bearing_differences.F)

## Summarise mean differences
bearings_between.summary.F <- summary(bearings_between.lm.F)

bearings_between.out.F <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between.lm.F"),
  intercept = bearings_between.summary.F$coefficients[1,1],
  intercept_sd = bearings_between.summary.F$coefficients[1,2],
  t(bearings_between.summary.F$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )


bearings_between.out.F %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))

## MALES ##
bearings_between.lm.M <- lm(bearing_diff ~ comp_ind, data = bearing_differences.M)

## Summarise mean differences
bearings_between.summary.M <- summary(bearings_between.lm.M)

bearings_between.out.M <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between.lm.M"),
  intercept = bearings_between.summary.M$coefficients[1,1],
  intercept_sd = bearings_between.summary.M$coefficients[1,2],
  t(bearings_between.summary.M$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )


bearings_between.out.M %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))



###### Between-individual, within seasons  -------------------------------

## FEMALES ##
bearings_between_withinseason.lmm.F <- lmer(bearing_diff ~ comp_ind + (1|season_1), 
                                          data = bearing_differences.filter.F)

## Extract mean differences from lmm
bearings_between_withinseason.summary.F <- summary(bearings_between_withinseason.lmm.F)

# Calculate individual repeatability
bearings_between_withinseason.varcor.F <- as.data.frame(VarCorr(bearings_between_withinseason.lmm.F))

ind_rpt.intercept.F <- bearings_between_withinseason.varcor.F$vcov[1] / (bearings_between_withinseason.varcor.F$vcov[1] + bearings_between_withinseason.varcor.F$vcov[2])


# Make summary output
bearings_between_withinseason.out.F <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between_withinseason.lmm.F"),
  intercept = bearings_between_withinseason.summary.F$coefficients[1,1],
  intercept_sd = bearings_between_withinseason.summary.F$coefficients[1,2],
  t(bearings_between_withinseason.summary.F$coefficients[2,]),
  intercept_var = bearings_between_withinseason.varcor.F$vcov[1],
  slope_var = NA,
  resid_var = bearings_between_withinseason.varcor.F$vcov[2]
) )


bearings_between_withinseason.out.F %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))


## MALES ##

bearings_between_withinseason.lmm.M <- lmer(bearing_diff ~ comp_ind + (1|season_1), 
                                            data = bearing_differences.filter.M)

## Extract mean differences from lmm
bearings_between_withinseason.summary.M <- summary(bearings_between_withinseason.lmm.M)

# Calculate individual repeatability
bearings_between_withinseason.varcor.M <- as.data.frame(VarCorr(bearings_between_withinseason.lmm.M))

ind_rpt.intercept.M <- bearings_between_withinseason.varcor.M$vcov[1] / (bearings_between_withinseason.varcor.M$vcov[1] + bearings_between_withinseason.varcor.M$vcov[2])


# Make summary output
bearings_between_withinseason.out.M <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "bearings_between_withinseason.lmm.M"),
  intercept = bearings_between_withinseason.summary.M$coefficients[1,1],
  intercept_sd = bearings_between_withinseason.summary.M$coefficients[1,2],
  t(bearings_between_withinseason.summary.M$coefficients[2,]),
  intercept_var = bearings_between_withinseason.varcor.M$vcov[1],
  slope_var = NA,
  resid_var = bearings_between_withinseason.varcor.M$vcov[2]
) )


bearings_between_withinseason.out.M %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))




# Plots to output -------------------------------------------------------------


### GPS & bearings from colony ----------------------------------------------

png(paste0("Figures/", my_species, "_", colony_exp, "GPS_and_bearings.png"),
    width = 14, height = 7, units = "in", res = 600)
ggarrange(gps_tracks.plot, bearings.plot,
                  ncol = 2)
dev.off()

### Comparison of bearing differences ---------------------------------------

bearings_within.plot2 <- bearings_within.plot +
  theme(plot.caption = element_blank())

bearings_between.plot2 <- bearings_between.plot +
  facet_wrap(.~comp_ind, ncol = 1, labeller = labeller(comp_ind = c(
    "within_ring" = "Within individual, between seasons",
    "between_ring" = "Between individual, between seasons" ))) +
  theme(plot.caption = element_blank(),
        axis.title.y = element_blank())

bearings_between_withinseason.plot2 <- bearings_between_withinseason.plot +
  facet_wrap(.~comp_ind, ncol = 1, labeller = labeller(comp_ind = c(
    "within_ring" = "Within individual, within season",
    "between_ring" = "Between individual, within season" ))) +
  theme(plot.caption = element_blank(),
        axis.title.y = element_blank())

png(paste0("Figures/", my_species, "_", colony_exp, "bearing_differences.png"),
    width = 15, height = 7, units = "in", res = 600)
ggarrange(bearings_within.plot2, bearings_between.plot2, bearings_between_withinseason.plot2,
          ncol = 3)
dev.off()

## Output statistical values -----------------------------------------------

# Bind outputs
bearings_output <- rbind(bearings_within.out, 
                         bearings_between.out,
                         bearings_between_withinseason.out)

# Add to file - replace any rows that already exist for this colony/species combo
bearings_output.existing <- read.csv("Data_outputs/bearings_model_summary.csv")
bearings_output.filtered <- bearings_output.existing %>%
  filter(!(species == my_species & colony == colony_exp))

bearings_output.save = rbind(bearings_output.filtered, bearings_output)

write.csv(bearings_output.save, "Data_outputs/bearings_model_summary.csv", row.names = F)


## With sex
bearings_output.sex <- rbind(bearings_between.out.F %>% mutate(sex = "F"),
                             bearings_between.out.M %>% mutate(sex = "M"),
                             bearings_between_withinseason.out.F %>% mutate(sex = "F"),
                             bearings_between_withinseason.out.M %>% mutate(sex = "M")) %>%
  relocate(sex, .after = colony)

# Add to file - replace any rows that already exist for this colony/species combo
bearings_output.sex.existing <- read.csv("Data_outputs/bearings_model_summary_bySex.csv")
bearings_output.sex.filtered <- bearings_output.sex.existing %>%
  filter(!(species == my_species & colony == colony_exp))

bearings_output.sex.save = rbind(bearings_output.sex.filtered, bearings_output.sex)

write.csv(bearings_output.sex.save, "Data_outputs/bearings_model_summary_bySex.csv", row.names = F)


### Repeatability values ----
rpt_df <- data.frame(colony = colony_exp, species = my_species, 
                     med_R.bearing = med_R.bearing, lci_R.bearing = ci_R.bearing[1],
                     uci_R.bearing = ci_R.bearing[2])
rownames(rpt_df) <- NULL

rpt_df.existing <- read.csv("Data_outputs/bearings_rpt.csv")
rpt_df.filtered <- rpt_df.existing %>%
  filter(!(species == my_species & colony == colony_exp))

rpt_df.save = rbind(rpt_df.filtered, rpt_df)

write.csv(rpt_df.save, "Data_outputs/bearings_rpt.csv", row.names = F)


# ______________________________ ####
# Approach 2: Comparing trips through Dynamic Time Warping (DTW) -------

rm(list = ls(pattern = "bearings"))

### Calculate DTW -------------------------------

# proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# 
# # Generate all possible trip combinations 
# trip_combinations <- combn(unique(my_gps$ID), 2, simplify = FALSE)
# trip_combinations.df <- data.frame(do.call("rbind", trip_combinations))
# colnames(trip_combinations.df) <- c("tripID_1", "tripID_2")
# 
# # Within-individual dataset
# dtw_within.df <- trip_combinations.df %>%
#   separate(tripID_1, into = "ring_1", sep = "_", remove = FALSE) %>%
#   separate(tripID_2, into = "ring_2", sep = "_", remove = FALSE) %>%
#   filter(ring_1 == ring_2) %>%
#   select(-c(ring_1, ring_2)) %>%
#   # Ensure comparisons are within same phase
#   left_join(my_meta %>% select(ID, phase) %>% 
#             rename(tripID_1 = ID),
#           by = "tripID_1") %>%
#   rename(phase_1 = phase) %>%
#   left_join(my_meta %>% select(ID, phase) %>% 
#               rename(tripID_2 = ID),
#             by = "tripID_2") %>%
#   rename(phase_2 = phase) %>%
#   mutate(phase_comp = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")) %>%
#   filter(phase_comp == "same_phase") %>%
#   select(-c(phase_comp, phase_1, phase_2))
# 
# ## Run through the comparisons
# dtw_within.list <- vector(mode = "list", length = nrow(dtw_within.df))
# 
# for (i in 1:nrow(dtw_within.df)) {
#   
#   print(paste0("Processing iteration ", i, " of ", nrow(dtw_within.df), "."))
#   
#   trip1 <- convert_to_spdf(subset(my_gps, ID == dtw_within.df$tripID_1[i]), proj.dec)
#   trip2 <- convert_to_spdf(subset(my_gps, ID == dtw_within.df$tripID_2[i]), proj.dec)
#   
#   trip1_mat <- as.matrix(trip1@coords)
#   trip2_mat <- as.matrix(trip2@coords)
#   
#   dtw_out <- dtw(trip1_mat, trip2_mat, dist.method = "Euclidean")$distance
#   
#   output <- data.frame(trip1 = dtw_within.df$tripID_1[i],
#                        trip2 = dtw_within.df$tripID_2[i],
#                        dtw = dtw_out)
# 
#   dtw_within.list[[i]] <- output 
#   
# }
# 
# dtw_within.output <- do.call("rbind", dtw_within.list)
# 
# dtw_within.output %<>% 
#   separate(trip1, into = c("ring_1", "season_1", "loc", "tripID_1"), sep = "_") %>%
#   separate(trip2, into = c("ring_2", "season_2", "loc", "tripID_2"), sep = "_") %>%
#   select(-loc) %>%
#   mutate(tripID_1 = paste0("trip_", tripID_1),
#          tripID_2 = paste0("trip_", tripID_2),
#          comp_ind = "within_ring",
#          comp_season = ifelse(season_1 == season_2, "within_season", "between_season"))
# 
# 
# # Between-individual dataset
# 
# dtw_between.df <- trip_combinations.df %>%
#   separate(tripID_1, into = "ring_1", sep = "_", remove = FALSE) %>%
#   separate(tripID_2, into = "ring_2", sep = "_", remove = FALSE) %>%
#   filter(ring_1 != ring_2) %>%
#   select(-c(ring_1, ring_2)) %>%
#   # Get phase and filter for same breeding phase
#   left_join(my_meta %>% select(ID, phase) %>% 
#               rename(tripID_1 = ID),
#             by = "tripID_1") %>%
#   rename(phase_1 = phase) %>%
#   left_join(my_meta %>% select(ID, phase) %>% 
#               rename(tripID_2 = ID),
#             by = "tripID_2") %>%
#   rename(phase_2 = phase) %>%
#   mutate(phase_comp = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")) %>%
#   filter(phase_comp == "same_phase") %>%
#   select(-c(phase_comp, phase_1, phase_2))
# 
# 
# ## Downsample the between-individual data
# set.seed(817)
# n = 10
# 
# dtw_between <- dtw_between.df %>% 
#   sample_n(size = n * nrow(dtw_within.df), replace = FALSE)
# 
# 
# ## Calculate the comparisons
# dtw_between.list <- vector(mode = "list", length = nrow(dtw_between))
# 
# for (i in 1:nrow(dtw_between)) {
#   
#   print(paste0("Processing iteration ", i, " of ", nrow(dtw_between), "."))
#   
#   trip1 <- convert_to_spdf(subset(my_gps, ID == dtw_between$tripID_1[i]), proj.dec)
#   trip2 <- convert_to_spdf(subset(my_gps, ID == dtw_between$tripID_2[i]), proj.dec)
#   
#   trip1_mat <- as.matrix(trip1@coords)
#   trip2_mat <- as.matrix(trip2@coords)
#   
#   dtw_out <- dtw(trip1_mat, trip2_mat, dist.method = "Euclidean")$distance
#   
#   output <- data.frame(trip1 = dtw_between$tripID_1[i],
#                        trip2 = dtw_between$tripID_2[i],
#                        dtw = dtw_out)
#   
#   dtw_between.list[[i]] <- output 
#   
# }
# 
# dtw_between.output <- do.call("rbind", dtw_between.list)
# 
# dtw_between.output %<>% 
#   separate(trip1, into = c("ring_1", "season_1", "loc", "tripID_1"), sep = "_") %>%
#   separate(trip2, into = c("ring_2", "season_2", "loc", "tripID_2"), sep = "_") %>%
#   select(-loc) %>%
#   mutate(tripID_1 = paste0("trip_", tripID_1),
#          tripID_2 = paste0("trip_", tripID_2),
#          comp_ind = "between_ring",
#          comp_season = ifelse(season_1 == season_2, "within_season", "between_season"))
# 
# 
# # Bind the two dataframes together
# dtw.df <- rbind(dtw_within.output, dtw_between.output)
# 
# save(dtw.df, file = paste0("Data_outputs/", my_species, "_", colony_exp, "_dtw_comparison.RData"))




### Prepare model data ------------------------------------------------------

load(file = paste0("Data_outputs/", my_species, "_", colony_exp, "_dtw_comparison.RData"))

## Merge sex data
load("Data_inputs/ind_meta.RData")

dtw.df.sex <- merge(dtw.df, all_ind.meta %>% select(-c(birthYr, species)) %>% rename(ring_1 = ring), 
                                 by = "ring_1", all.x = T)

#### Within individuals, across seasons -------------------------------

dtw_within.loglmm <- lmer(log(dtw) ~ comp_season + (1+comp_season|ring_1), data = dtw.df %>% filter(comp_ind == "within_ring"))

## Extract mean differences and confint from log LMM
dtw_within.summary <- summary(dtw_within.loglmm)

ci_between <- exp(fixef(dtw_within.loglmm)["(Intercept)"] + c(-1, 1) * 1.96 * summary(dtw_within.loglmm)$coefficients["(Intercept)", "Std. Error"])
ci_within <- exp(fixef(dtw_within.loglmm)["(Intercept)"] + fixef(dtw_within.loglmm)["comp_seasonwithin_season"] + c(-1, 1) * 1.96 * summary(dtw_within.loglmm)$coefficients["comp_seasonwithin_season", "Std. Error"])

dtw_within.means <- data.frame(comp_season = c("between_season", "within_season"),
                               mean = exp(c(dtw_within.summary$coefficients[1,1], dtw_within.summary$coefficients[1,1] + dtw_within.summary$coefficients[2,1])),
                               ci_low = c(ci_between[1], ci_within[1]),
                               ci_high = c(ci_between[2], ci_within[2]))

# Plot #
dtw_within_ind.plot <- ggplot() +
  geom_histogram(data = dtw.df %>% filter(comp_ind == "within_ring"), aes(x = dtw)) +
  geom_segment(data = dtw_within.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = dtw_within.means, aes(x = ci_low, xend = ci_low, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = dtw_within.means, aes(x = ci_high, xend = ci_high, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_season, ncol = 1, labeller = labeller(comp_season = c(
    "within_season" = "Within seasons",
    "between_season" = "Between seasons" ))) +
  labs(y = "Frequency", x = "DTW distance (track dissimilarity)",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear mixed effects model with logged response.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0),
        strip.background = element_rect(fill ="grey90"),
        strip.text = element_text(face = "bold"))


# Calculate individual repeatability
dtw_within.varcor <- as.data.frame(VarCorr(dtw_within.loglmm))

ind_rpt.intercept <- dtw_within.varcor$vcov[1] / (dtw_within.varcor$vcov[1] + dtw_within.varcor$vcov[4])
ind_rpt.slope <- dtw_within.varcor$vcov[2] / (dtw_within.varcor$vcov[2] + dtw_within.varcor$vcov[4])


# Make a summary output
dtw_within.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_within.lmm"),
  intercept = dtw_within.summary$coefficients[1,1],
  intercept_sd = dtw_within.summary$coefficients[1,2],
  t(dtw_within.summary$coefficients[2,]),
  intercept_var = dtw_within.varcor$vcov[1],
  slope_var = dtw_within.varcor$vcov[2],
  resid_var = dtw_within.varcor$vcov[4]
) )


dtw_within.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))


#### Between individuals, across seasons -----------------------------------------------------

dtw_between.loglm <- lm(log(dtw) ~ comp_ind, data = dtw.df)

dtw_between.summary <- summary(dtw_between.loglm)

ci_log <- confint(dtw_between.loglm)
ci_between <- exp(ci_log["(Intercept)",])
ci_within <- exp(ci_log["(Intercept)", ] + ci_log["comp_indwithin_ring", ])

dtw_between.means <- data.frame(comp_season = c("between_season", "within_season"),
                                mean = c(exp(coef(dtw_between.loglm)["(Intercept)"]),
                                         exp(coef(dtw_between.loglm)["(Intercept)"] +    coef(dtw_between.loglm)["comp_indwithin_ring"])),
                                ci_low = c(ci_between[1], ci_within[1]),
                                ci_high = c(ci_between[2], ci_within[2]))

# Plot #
dtw_between.plot <- ggplot() +
  geom_histogram(data = dtw.df %>% filter(comp_ind == "within_ring"), aes(x = dtw)) +
  geom_segment(data = dtw_within.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = dtw_within.means, aes(x = ci_low, xend = ci_low, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = dtw_within.means, aes(x = ci_high, xend = ci_high, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_season, ncol = 1, labeller = labeller(comp_season = c(
    "within_season" = "Within Season",
    "between_season" = "Between Season" ))) +
  labs(y = "Frequency", x = "DTW distance (track dissimilarity)",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear model with logged response.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0),
        strip.background = element_rect(fill ="grey90"),
        strip.text = element_text(face = "bold"))



# Make summary output
dtw_between.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between.lm"),
  intercept = dtw_between.summary$coefficients[1,1],
  intercept_sd = dtw_between.summary$coefficients[1,2],
  t(dtw_between.summary$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )

dtw_between.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))


#### Between individual, Within seasons -------

# Need to ensure that only calculate when have multiple birds per season
dtw.df.filter <- dtw.df %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")

dtw_between_withinseason.loglmm <- lmer(log(dtw) ~ comp_ind + (1|season_1), data = dtw.df.filter)

## Extract mean differences and confint from log LMM
dtw_between_withinseason.summary <- summary(dtw_between_withinseason.loglmm)

ci_between <- exp(fixef(dtw_between_withinseason.loglmm)["(Intercept)"] + c(-1, 1) * 1.96 * summary(dtw_between_withinseason.loglmm)$coefficients["(Intercept)", "Std. Error"])
ci_within <- exp(fixef(dtw_between_withinseason.loglmm)["(Intercept)"] + fixef(dtw_between_withinseason.loglmm)["comp_indwithin_ring"] + c(-1, 1) * 1.96 * summary(dtw_between_withinseason.loglmm)$coefficients["comp_indwithin_ring", "Std. Error"])

dtw_between_withinseason.means <- data.frame(comp_ind = c("between_ring", "within_ring"),
                                             mean = exp(c(dtw_between_withinseason.summary$coefficients[1,1], dtw_between_withinseason.summary$coefficients[1,1] + dtw_between_withinseason.summary$coefficients[2,1])),
                                             ci_low = c(ci_between[1], ci_within[1]),
                                             ci_high = c(ci_between[2], ci_within[2]))

# Plot #
dtw_between_withinseason.plot <- ggplot() +
  geom_histogram(data = dtw.df %>% filter(comp_season == "within_season"), aes(x = dtw)) +
  geom_segment(data = dtw_between_withinseason.means, aes(x = mean, xend = mean, y = 0, yend = Inf), 
               color = "red", linewidth = 1) +
  geom_segment(data = dtw_between_withinseason.means, aes(x = ci_low, xend = ci_low, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = dtw_between_withinseason.means, aes(x = ci_high, xend = ci_high, y = 0, yend = Inf), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(.~comp_ind, ncol = 1, labeller = labeller(comp_ind = c(
    "between_ring" = "Between individuals",
    "within_ring" = "Within individuals" ))) +
  labs(y = "Frequency", x = "DTW distance (track dissimilarity)",
       caption = "Red lines indicate mean (solid) \u00B1 standard error (dashed) taken from linear mixed effects model with logged response.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0),
        strip.background = element_rect(fill ="grey90"),
        strip.text = element_text(face = "bold"))



# Calculate individual repeatability
dtw_between_withinseason.varcor <- as.data.frame(VarCorr(dtw_between_withinseason.loglmm))

ind_rpt.intercept <- dtw_between_withinseason.varcor$vcov[1] / (dtw_between_withinseason.varcor$vcov[1] + dtw_between_withinseason.varcor$vcov[2])


# Make summary output
dtw_between_withinseason.out <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between_withinseason.lmm"),
  intercept = dtw_between_withinseason.summary$coefficients[1,1],
  intercept_sd = dtw_between_withinseason.summary$coefficients[1,2],
  t(dtw_between_withinseason.summary$coefficients[2,]),
  intercept_var = dtw_between_withinseason.varcor$vcov[1],
  slope_var = NA,
  resid_var = dtw_between_withinseason.varcor$vcov[2]
) )


dtw_between_withinseason.out %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))


### By sex -----------------------

dtw.df.F <- dtw.df.sex %>% filter(sex == "F")
dtw.df.M <- dtw.df.sex %>% filter(sex == "M")


# Need to ensure that only calculate when have multiple birds per season
dtw.df.filter.F <- dtw.df.F %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")

dtw.df.filter.M <- dtw.df.M %>% group_by(season_1) %>% mutate(n_birds = n_distinct(ring_1)) %>%
  filter(n_birds > 1) %>% select(-n_birds) %>% filter(comp_season == "within_season")


##### Between individuals, across seasons -----------------------------------------------------

## FEMALE ##
dtw_between.loglm.F <- lm(log(dtw) ~ comp_ind, data = dtw.df.F)
dtw_between.summary.F <- summary(dtw_between.loglm.F)

# Make summary output
dtw_between.out.F <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between.lm.F"),
  intercept = dtw_between.summary.F$coefficients[1,1],
  intercept_sd = dtw_between.summary.F$coefficients[1,2],
  t(dtw_between.summary.F$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )

dtw_between.out.F %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))


## MALE ##
dtw_between.loglm.M <- lm(log(dtw) ~ comp_ind, data = dtw.df.M)
dtw_between.summary.M <- summary(dtw_between.loglm.M)

# Make summary output
dtw_between.out.M <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between.lm.M"),
  intercept = dtw_between.summary.M$coefficients[1,1],
  intercept_sd = dtw_between.summary.M$coefficients[1,2],
  t(dtw_between.summary.M$coefficients[2,]),
  intercept_var = NA,
  slope_var = NA,
  resid_var = NA
) )

dtw_between.out.M %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(df = NA) %>%
  relocate(df, .before = t_value) %>%
  mutate(across(intercept:resid_var, as.numeric))



##### Between individual, Within seasons -------

## FEMALE ##
dtw_between_withinseason.loglmm.F <- lmer(log(dtw) ~ comp_ind + (1|season_1), data = dtw.df.F %>% filter(comp_season == "within_season"))

## Extract mean differences and confint from log LMM
dtw_between_withinseason.summary.F <- summary(dtw_between_withinseason.loglmm.F)

# Calculate individual repeatability
dtw_between_withinseason.varcor.F <- as.data.frame(VarCorr(dtw_between_withinseason.loglmm.F))

ind_rpt.intercept.F <- dtw_between_withinseason.varcor.F$vcov[1] / (dtw_between_withinseason.varcor.F$vcov[1] + dtw_between_withinseason.varcor.F$vcov[2])

# Make summary output
dtw_between_withinseason.out.F <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between_withinseason.lmm.F"),
  intercept = dtw_between_withinseason.summary.F$coefficients[1,1],
  intercept_sd = dtw_between_withinseason.summary.F$coefficients[1,2],
  t(dtw_between_withinseason.summary.F$coefficients[2,]),
  intercept_var = dtw_between_withinseason.varcor.F$vcov[1],
  slope_var = NA,
  resid_var = dtw_between_withinseason.varcor.F$vcov[2]
) )


dtw_between_withinseason.out.F %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))



## FEMALE ##
dtw_between_withinseason.loglmm.M <- lmer(log(dtw) ~ comp_ind + (1|season_1), data = dtw.df.M %>% filter(comp_season == "within_season"))

## Extract mean differences and confint from log LMM
dtw_between_withinseason.summary.M <- summary(dtw_between_withinseason.loglmm.M)

# Calculate individual repeatability
dtw_between_withinseason.varcor.M <- as.data.frame(VarCorr(dtw_between_withinseason.loglmm.M))

ind_rpt.intercept.M <- dtw_between_withinseason.varcor.M$vcov[1] / (dtw_between_withinseason.varcor.M$vcov[1] + dtw_between_withinseason.varcor.M$vcov[2])

# Make summary output
dtw_between_withinseason.out.M <- data.frame(cbind(
  cbind(species = my_species,
        colony = colony_exp,
        model = "dtw_between_withinseason.lmm.M"),
  intercept = dtw_between_withinseason.summary.M$coefficients[1,1],
  intercept_sd = dtw_between_withinseason.summary.M$coefficients[1,2],
  t(dtw_between_withinseason.summary.M$coefficients[2,]),
  intercept_var = dtw_between_withinseason.varcor.M$vcov[1],
  slope_var = NA,
  resid_var = dtw_between_withinseason.varcor.M$vcov[2]
) )


dtw_between_withinseason.out.M %<>% 
  rename(estimate = Estimate, se = Std..Error, t_value = t.value, p_value = Pr...t..) %>%
  mutate(across(intercept:resid_var, as.numeric))








# Output statistical values -----------------------------------------------

# Bind outputs
dtw_output <- rbind(dtw_within.out, 
                    dtw_between.out,
                    dtw_between_withinseason.out)

# Add to file - replace any rows that already exist for this colony/species combo
dtw_output.existing <- read.csv("Data_outputs/dtw_model_summary.csv")
dtw_output.filtered <- dtw_output.existing %>%
  filter(!(species == my_species & colony == colony_exp))

dtw_output.save = rbind(dtw_output.filtered, dtw_output)

write.csv(dtw_output.save, "Data_outputs/dtw_model_summary.csv", row.names = F)



# With sex
dtw_output.sex <- rbind(dtw_between.out.F %>% mutate(sex = "F"),
                        dtw_between.out.M %>% mutate(sex = "M"),
                        dtw_between_withinseason.out.F %>% mutate(sex = "F"),
                        dtw_between_withinseason.out.M %>% mutate(sex = "M")) %>%
  relocate(sex, .after = colony)

# Add to file - replace any rows that already exist for this colony/species combo
dtw_output.sex.existing <- read.csv("Data_outputs/dtw_model_summary_bySex.csv")
dtw_output.sex.filtered <- dtw_output.sex.existing %>%
  filter(!(species == my_species & colony == colony_exp))

dtw_output.sex.save = rbind(dtw_output.sex.filtered, dtw_output.sex)

write.csv(dtw_output.sex.save, "Data_outputs/dtw_model_summary_bySex.csv", row.names = F)



# ______________________________ ####
# Approach 3: Similarity in patch use -------

### Calculate bearings to foraging patches ----------------------------------

load(paste0("Data_outputs/", my_species, "_", colony_exp, "_labelledHMM.RData"))

my_gps.forage <- my_gps %>%
  group_by(ring, ID) %>%
  mutate(state_chg = State != lag(State, default = first(State))) %>%  # Identify changes in state
  mutate(state_grp = cumsum(state_chg)) %>% 
  group_by(State) %>%
  mutate(State_ID = paste(State, match(state_grp, unique(state_grp)), sep = "_")) %>%
  ungroup() %>%
  select(-state_chg, -state_grp) 

forage_bouts <- my_gps.forage %>%
  group_by(ring, ID, State_ID) %>%
  summarise(stt = datetime[1],
         season = season[1],
         med_lon = median(longitude, na.rm = T),
         med_lat = median(latitude, na.rm = T)) %>%
  filter(grepl("Search_1", State_ID)) %>%
  mutate(
    bearing_col = pmap_dbl(
      list(med_lon, med_lat),
      ~ bearing(c(col_coords$longitude, col_coords$latitude), c(..1, ..2)) ))

## Convert bearings to radians
forage_bouts %<>% mutate(bearing_rad = deg2rad(bearing_col))


### Repeatability of foraging patch bearings ------------------------------------

load("Data_inputs/ind_meta.RData")

# Integrate meta
forage_bearings.rpt <- merge(forage_bouts, all_ind.meta %>% filter(species == toupper(my_species)), by = "ring")

# Fit the model
forage_bearings_rpt.brms <- brm(
  formula = bf(bearing_rad ~ (1|ring), family = von_mises()),
  data = forage_bearings.rpt,
  chains = 4, iter = 2000, cores = 4)

summary(forage_bearings_rpt.brms)

posterior_samples <- posterior_samples(forage_bearings_rpt.brms, pars = c("sd_ring__Intercept", "kappa"))

sd_individual <- posterior_samples$sd_ring__Intercept  
kappa <- posterior_samples$kappa  
var_within <- 1 / kappa  

repeatability <- sd_individual^2 / (sd_individual^2 + var_within)


#### Output bearing repeatability summary statistics -----
med_R.bearing <- median(repeatability) 
ci_R.bearing <- quantile(repeatability, c(0.025, 0.975))  

rpt_df <- data.frame(colony = colony_exp, species = my_species, 
                     med_R.bearing = med_R.bearing, lci_R.bearing = ci_R.bearing[1],
                     uci_R.bearing = ci_R.bearing[2])
rownames(rpt_df) <- NULL

rpt_df.existing <- read.csv("Data_outputs/forage_bearings_rpt.csv")
rpt_df.filtered <- rpt_df.existing %>%
  filter(!(species == my_species & colony == colony_exp))

rpt_df.save = rbind(rpt_df.filtered, rpt_df)

write.csv(rpt_df.save, "Data_outputs/forage_bearings_rpt.csv", row.names = F)




# ______________________________ ####
# Plot GPS trips by season ----------------------------------------------------------

# # Get colony locations
# if(colony == "cro" | colony == "ker") { colony.df = data.frame(Lon = 51.706972, Lat = -46.358639) }
# if(colony == "bi") { colony.df = data.frame(Lon = -38.05, Lat = -54) }
# 
# # Set projections
# proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# proj.utm <- paste0("+proj=laea +lat_0=", colony.df$Lat[1], " +lon_0=", colony.df$Lon[1],
#                    " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# ## Colony projection (for plotting)
# coordinates(colony.df) <- ~Lon+Lat
# proj4string(colony.df) <- proj.dec
# colony.sp <- spTransform(colony.df, CRS(proj.utm))
# colony.df <- as.data.frame(colony.sp)
# 
# ## Get the base world plot
# world <- ne_countries(scale = "medium", returnclass = "sf")
# world2 <- sf::st_transform(world, crs = proj.utm)
# 
# ## Set projection of GPS data
# test_trips.proj <- test_trips
# coordinates(test_trips.proj) <- ~ longitude+latitude
# proj4string(test_trips.proj) <- proj.dec
# 
# test_trips.proj.sp <- spTransform(test_trips.proj, CRS(proj.utm))
# test_trips.proj.df <- as.data.frame(test_trips.proj.sp)
# 
# 
# padding <- (max(test_trips.proj.df$coords.x1) - min(test_trips.proj.df$coords.x1)) * 0.2 
#   
# my_xmin <- min(test_trips.proj.df$coords.x1) - padding
# my_xmax <- max(test_trips.proj.df$coords.x1) + padding
# my_ymin <- min(test_trips.proj.df$coords.x2) - padding
# my_ymax <- max(test_trips.proj.df$coords.x2) + padding
#   
# ggplot(data = world2) + 
#     geom_sf(fill = "cadetblue", colour = "grey") +
#     coord_sf(crs = proj.utm, xlim = c(my_xmin, my_xmax), ylim = c(my_ymin, my_ymax), 
#              label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
#     geom_path(aes(x = coords.x1, y = coords.x2, colour = boutID), size = 1, dat = test_trips.proj.df) +
#     scale_colour_viridis_d(name = NULL) +
#     labs(title = paste0("Ring: ", test_trips.proj.df$ring[1])) +
#     facet_grid(.~season) + 
#     theme(panel.background = element_rect(fill = "white"), 
#           panel.grid.major = element_line(colour = "grey80"),
#           panel.border = element_rect(colour = "black", fill = NA),
#           axis.text.x.top = element_blank(), 
#           axis.title.x.top = element_blank(),
#           axis.ticks.x.top = element_blank(),
#           axis.title.x.bottom = element_blank(),
#           axis.text.y.right = element_blank(), 
#           axis.title.y.right = element_blank(),
#           axis.title.y.left = element_blank(),
#           plot.title = element_text(size = 12, face = "italic"),
#           legend.position = "none", 
#           strip.text = element_text(face = "bold")) 
  