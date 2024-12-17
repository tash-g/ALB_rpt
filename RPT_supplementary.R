
## ---------------------------
##
## Script name: RPT_supplementary
##
## Purpose of script: Supplementary analyses for repeatability paper
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-12-09
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("ALB_FOR_functions.R")
source("RPT_functions.R")

# Define the packages
packages <- c("dplyr", "magrittr", "ggplot2", "lme4", "rptR", "gridExtra", "tidyr")

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


# Compare repeatability between phases and colonies -----------------------

load("Data_inputs/repeatability_data.RData")
rep_bba <- rep_data %>% filter(species == "bba")
rep_waal <- rep_data %>% filter(species == "waal")

rep_bba %<>% mutate(travRestRatio = transit_time.hrs/rest_time.hrs)
rep_waal %<>% mutate(travRestRatio = transit_time.hrs/rest_time.hrs)

# Trying different response variables and examining repeatability
cov <- c("duration.hrs", "transitions", "transit_time.pct", "rest_time.pct", 
         "total_distance.km", "max_distance.km", "travRestRatio")

source("RPT_functions.R")


rpt.grouped_list <- vector(mode = "list", length = length(cov))

for (i in 1:length(cov)) {
  
  print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M"), " Processing covariate ", i, " of ", length(cov), ": ", cov[i]))
  
  formula <- as.formula(paste(cov[i], "~ (1|ring) + (1|season) + (1|ring:season) + age"))
  
  # BBAL 
  
  ## Breeding phase
  bba_phase.rpt <- calc_rpt_grouped(rep_bba, formula, "phase", c("incubation", "brooding"), "BBAL") 
  
  ## Colony
  bba_col.rpt <- calc_rpt_grouped(rep_bba, formula, "colony", c("birdis", "kerguelen"), "BBAL") 
  
  bba_rpt.grouped <- cbind( bba_phase.rpt, 
                            bba_col.rpt %>% select(-c(species, group)) %>% rename(R_french = R_kerguelen))
  
  
  # WAAL 
  
  ## Breeding phase
  waal_phase.rpt <- calc_rpt_grouped(rep_waal, formula, "phase", c("incubation", "brooding"), "WAAL") 
  
  ## Colony
  waal_col.rpt <- calc_rpt_grouped(rep_waal, formula, "colony", c("birdis", "crozet"), "WAAL") 
  
  waal_rpt.grouped <- cbind( waal_phase.rpt, 
                             waal_col.rpt %>% select(-c(species, group)) %>% rename(R_french = R_crozet))
  
  # Output
  rpt.grouped <- rbind(bba_rpt.grouped, waal_rpt.grouped)
  rpt.grouped %<>%
    mutate(cov = cov[i]) %>%
    relocate(cov, .before = species)
  
  rpt.grouped %<>% 
    pivot_wider(
      names_from = species,             
      values_from = c(R_incubation, R_brooding, R_birdis, R_french),      
      names_glue = "{species}_{.value}" 
    ) %>% 
    rename(parameter = cov) %>%
    relocate(parameter, group, BBAL_R_incubation, BBAL_R_brooding, 
             BBAL_R_birdis, BBAL_R_french)
  
  rpt.grouped_list[[i]] <- rpt.grouped
  
}

rpt.grouped_df <- do.call("rbind", rpt.grouped_list)
readr::write_excel_csv(rpt.grouped_df, "Data_outputs/SUPPLEMENTARY_grouped_repeatability.csv")


## Sample sizes
rpt_ind.data <- rep_data %>% group_by(ring) %>% select(species, ring, colony, phase) %>% distinct()
table(rpt_ind.data$colony, rpt_ind.data$species)
table(rpt_ind.data$phase, rpt_ind.data$species)

rpt_ind.data.gps <- rep_data %>% filter(!is.na(max_distance.km)) %>% group_by(ring) %>% mutate(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1) %>% ungroup() %>% select(-n_seasons) 
rpt_ind.data.gps %<>% group_by(ring) %>% select(ring, species, colony, phase) %>% distinct()

table(rpt_ind.data.gps$phase, rpt_ind.data.gps$species)
table(rpt_ind.data.gps$colony, rpt_ind.data.gps$species)

