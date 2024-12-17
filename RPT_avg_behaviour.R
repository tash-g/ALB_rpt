## ---------------------------
##
## Script name: RPT_avg_behaviour.R
##
## Purpose of script: Quantify repeatability in individual foraging behaviour
##
## Author: Dr. Natasha Gillies
##
## Created: 2024-08-05
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


# Load & process data ==========================================================

ind_files <- list.files("Data_inputs/", pattern = "individual_trips")
ind_list <- vector(mode = "list", length = 4)

for (i in 1:length(ind_files)) {
  
  my_data <- loadRData(paste0("Data_inputs/",ind_files[i]))
  ind_list[[i]] <- my_data
  
}

trip_data <- do.call("rbind", ind_list)

### Load metadata --------------------------------------------------------------

meta_files <- list.files("Data_inputs/", pattern = "individual_meta")
meta_data.A <- loadRData(paste0("Data_inputs/", meta_files[1]))
meta_data.B <- loadRData(paste0("Data_inputs/", meta_files[2])) %>% select(-colony) %>% 
  rename(birthYr = birth_yr) %>% relocate(ring, sex, species, birthYr)
meta_data <- rbind(meta_data.A, meta_data.B) %>% mutate(species = tolower(species))

## Bind it in
trip_data <- merge(trip_data, meta_data, by = c("ring", "species")) %>%
  mutate(birthYr = as.numeric(birthYr),
         age = season - birthYr) 

trip_data %<>% mutate(age = ifelse(age < 0, NA, age))
trip_data.full <- trip_data

# ______________________________ ####
# Get sample sizes -----------------------------------------------------------

### All birds -----------
trip_data %<>% filter(!is.na(age))

# Birds per colony/species/phase
ind_data <- trip_data %>% group_by(ring) %>% select(ring, species, colony, phase) %>% distinct()
table(ind_data$colony, ind_data$phase, ind_data$species)
table(ind_data$colony, ind_data$species)

## More than one season
ind_data.rep <- trip_data %>% group_by(ring) %>% mutate(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1) %>% ungroup() %>% select(-n_seasons) %>% select(-sex)
ind_data.rep %<>% group_by(ring) %>% select(ring, species, colony, phase) %>% distinct()

table(ind_data.rep$colony, ind_data.rep$phase, ind_data.rep$species)

#### GPS only -----
ind_data <- trip_data %>% filter(!is.na(max_distance.km)) %>% group_by(ring) %>% select(ring, species, colony, phase) %>% distinct()
table(ind_data$colony, ind_data$phase, ind_data$species)

## More than one season
ind_data.rep <- trip_data %>% filter(!is.na(max_distance.km)) %>% group_by(ring) %>% mutate(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1) %>% ungroup() %>% select(-n_seasons) %>% select(-sex)
ind_data.rep %<>% group_by(ring) %>% select(ring, species, colony, phase) %>% distinct()

table(ind_data.rep$colony, ind_data.rep$phase, ind_data.rep$species)

## Trips -----------

# Trips per colony/species/phase
table(trip_data$colony, trip_data$phase, trip_data$species)
table(trip_data$colony, trip_data$species)

trip_data2 <- trip_data %>% filter(!is.na(max_distance.km))
table(trip_data2$colony, trip_data2$phase, trip_data2$species)


# Trips per bird
trip_data %>% group_by(ring, colony, species) %>% summarise(n_trips = n()) %>%
  ungroup() %>%
  group_by(colony, species) %>% summarise(n_trips = mean(n_trips))

# Seasons per bird
trip_data %>% group_by(ring, colony, species) %>% summarise(n_seasons = n_distinct(season)) %>%
  ungroup() %>%
  group_by(colony, species) %>% summarise(n_seasons = mean(n_seasons))

# Seasons per colony
trip_data %>% group_by(colony, species) %>% summarise(n_seasons = n_distinct(season)) 


#### Birds with > 1 season -----------
# Isolate birds with > 1 season
rep_data <- trip_data %>% group_by(ring) %>% mutate(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1) %>% ungroup() %>% select(-n_seasons) %>% select(-sex)

save(rep_data, file = "Data_inputs/repeatability_data.RData")

# Birds per colony/species
ind_data2 <- rep_data %>% group_by(ring) %>% select(ring, species, colony) %>% distinct()
table(ind_data2$colony, ind_data2$species)

# Trips per colony/species
table(rep_data$colony, rep_data$species)

# Trips per bird
rep_data %>% group_by(ring, colony, species) %>% summarise(n_trips = n()) %>%
  ungroup() %>%
  group_by(colony, species) %>% summarise(n_trips = mean(n_trips))

# Seasons per bird
rep_data %>% group_by(ring, colony, species) %>% summarise(n_seasons = n_distinct(season)) %>%
  ungroup() %>%
  group_by(colony, species) %>% summarise(n_seasons = mean(n_seasons))

# Seasons per colony
trip_data %>% group_by(colony, species) %>% summarise(n_seasons = n_distinct(season)) 

# Number of between-year repeats per individual
trip_data %>%
  group_by(colony, species, ring) %>%
  summarise(n_seasons = n_distinct(season)) %>%
  ungroup() %>%
  mutate(rpt_category = case_when(
    n_seasons == 1 ~ "1 seasons",
    n_seasons == 2 ~ "2 seasons",
    n_seasons == 3 ~ "3 seasons",
    n_seasons > 3 ~ "4+ seasons"
  )) %>%
  group_by(colony, species, rpt_category) %>%
  summarise(num_birds = n()) %>%
  arrange(colony, species, rpt_category) %>%
  pivot_wider(
    names_from = rpt_category,       
    values_from = num_birds, 
    values_fill = 0)



# ______________________________ ####
# Examining repeatability ======================================================

load("Data_inputs/repeatability_data.RData")
rep_bba <- rep_data %>% filter(species == "bba")
rep_waal <- rep_data %>% filter(species == "waal")

rep_bba %<>% mutate(travRestRatio = transit_time.hrs/rest_time.hrs)
rep_waal %<>% mutate(travRestRatio = transit_time.hrs/rest_time.hrs)

# Trying different response variables and examining repeatability
cov <- c("duration.hrs", "transitions", "transit_time.pct", "rest_time.pct", 
         "total_distance.km", "max_distance.km", "travRestRatio")

rpt_output.list <- vector(mode = "list", length = length(cov))

for (i in 1:length(cov)) {
  
  print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M"), " Processing covariate ", i, " of ", length(cov), ": ", cov[i]))
  
  formula <- as.formula(paste(cov[i], "~ (1|ring) + (1|season) + (1|ring:season) + phase + age"))
  
  tok <- Sys.time() 
  
  repeatability_output.bba <- boot_rpt(1000, rep_bba.cov, formula)
  repeatability_output.waal <- boot_rpt(1000, rep_waal.cov, formula)
  
  tik <- Sys.time()
  print(paste0("Minutes to bootstrap repeatability: ", tik-tok))
  
  # Permutation test of repeatability
  tok <- Sys.time() 
  
  perm_output.bba <- perm_repeatability(1000, rep_bba.cov, formula, 
                                        repeatability_output.bba$est[1], 
                                        repeatability_output.bba$est[2],
                                        repeatability_output.bba$est[3])
  
  perm_output.waal <- perm_repeatability(1000, rep_waal.cov, formula, 
                                         repeatability_output.waal$est[1], 
                                         repeatability_output.waal$est[2],
                                         repeatability_output.waal$est[3])
  
  tik <- Sys.time()
  print(paste0("Minutes to complete permutation: ", tik-tok))
  
  ### Summarise repeatability output ----
  
  # BBAL
  
  ## Get variances
  rpt.bba <- lmer(formula, rep_bba.cov)
  vars.bba <- extract_var(rpt.bba)
  
  ## Summary table
  rpt_summary.bba <- cbind(vars.bba, repeatability_output.bba, perm_output.bba)
  rpt_summary.bba$group <- rownames(rpt_summary.bba)
  rpt_summary.bba$group = gsub("_var", "", rpt_summary.bba$group)
  rownames(rpt_summary.bba) <- NULL
  
  rpt_summary.bba %<>% 
    mutate(species = "BBAL",
           parameter = cov[i],
           R = paste0(signif(est, 2), " ± ", signif(se, 2)),
           group = gsub("se_", "", group)) %>%
    relocate(c(group, R), .before = se) %>%
    select(-c(est, se)) %>%
    pivot_wider(
      names_from = group,             
      values_from = c(R, ci, var, p),      
      names_glue = "{group}_{.value}" 
    ) %>% 
    relocate(species, parameter, annual_var, annual_R, annual_ci, annual_p,
             within_var, within_R, within_ci, within_p, 
             between_var, between_R, between_ci, between_p) %>%
    mutate(sample_size = paste0(nrow(rep_bba.cov), ", ", n_distinct(rep_bba.cov$ring)))
  
  # WAAL
  
  ## Get variances
  rpt.waal <- lmer(formula, rep_waal.cov)
  vars.waal <- extract_var(rpt.waal)
  
  ## Summary table
  rpt_summary.waal <- cbind(vars.waal, repeatability_output.waal, perm_output.waal)
  rpt_summary.waal$group <- rownames(rpt_summary.waal)
  rpt_summary.waal$group = gsub("_var", "", rpt_summary.waal$group)
  rownames(rpt_summary.waal) <- NULL
  
  rpt_summary.waal %<>% 
    mutate(species = "WAAL",
           parameter = cov[i],
           R = paste0(signif(est, 2), " ± ", signif(se, 2)),
           group = gsub("se_", "", group)) %>%
    relocate(c(group, R), .before = se) %>%
    select(-c(est, se)) %>%
    pivot_wider(
      names_from = group,             
      values_from = c(R, ci, var, p),      
      names_glue = "{group}_{.value}" 
    ) %>% 
    relocate(species, parameter, annual_var, annual_R, annual_ci, annual_p,
             within_var, within_R, within_ci, within_p, 
             between_var, between_R, between_ci, between_p) %>%
    mutate(sample_size = paste0(nrow(rep_waal.cov), ", ", n_distinct(rep_waal.cov$ring)))
  
  # Output summary
  rpt_output.list[[i]] <- rbind(rpt_summary.bba, rpt_summary.waal)
   
}


rpt_output.full <- do.call("rbind", rpt_output.list)
rpt_output.full %<>% arrange(species)


readr::write_excel_csv(rpt_output.full, "Data_outputs/repeatability_output.csv")




# ______________________________ ####
# Integrating foraging measures - PCA ==========================================

### Process the data --------
load("Data_inputs/repeatability_data.RData")

bbal.df <- rep_data %>% filter(species == "bba")
waal.df <- rep_data %>% filter(species == "waal")

# Packages
library(corrr); library(ggcorrplot); library(FactoMineR); library(factoextra)

# Check for missing values - these can bias PCA
colSums(is.na(bbal.df))
colSums(is.na(waal.df))

# Normalise the data
bbal_pca.df <- bbal.df %>%
  select(duration.hrs, transitions, transit_time.pct, rest_time.pct) %>%
  mutate(across(everything(), ~scale(.) %>% as.vector()))

waal_pca.df <- waal.df %>%
  select(duration.hrs, transitions, transit_time.pct, rest_time.pct) %>%
  mutate(across(everything(), ~scale(.) %>% as.vector()))

### Conduct the PCA ----

bbal.pca <- princomp(bbal_pca.df)
summary(bbal.pca)
loadings(bbal.pca)

waal.pca <- princomp(waal_pca.df)
summary(waal.pca)

# Look at importance of each PC
grid.arrange(
  fviz_eig(bbal.pca, addlabels = TRUE, main = "BBAL"),
  fviz_eig(waal.pca, addlabels = TRUE, main = "WAAL"),
  ncol = 2 )

# Biplot of attributes
grid.arrange(
  fviz_pca_var(bbal.pca, repel = TRUE, title = "BBAL"),
  fviz_pca_var(waal.pca, repel = TRUE, title = "WAAL"),
  ncol = 2 )

# Contribution of each variable to each component
fviz_cos2(bbal.pca, choice = "var", axes = 1:2)
fviz_cos2(waal.pca, choice = "var", axes = 1:2)

# Combined biplot and cosplot
fviz_pca_var(bbal.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)


### Extract PCA components -----
bbal_pca.df <- cbind(bbal.df %>% select(c(ring, species, colony, season, phase, sex, age)),
                     bbal.pca$scores)

waal_pca.df <- cbind(waal.df %>% select(c(ring, species, colony, season, phase, sex, age)),
                     waal.pca$scores)

library(brms)

rpt_brms.bba <- brm(
  mvbind(Comp.1, Comp.2) ~ sex + age + (1 | ring) + (1 | season) + (1|ring:season),
  data = bbal_pca.df,  
  family = gaussian(), 
  chains = 4,          
  cores = 4,           
  iter = 2000          
)

save(rpt_brms.bba, file = "Data_outputs/multivariate_rpt_model_bbal.RData")

rpt_brms.waal <- brm(
  mvbind(Comp.1, Comp.2) ~ sex + age + (1 | ring) + (1 | season) + (1|ring:season),
  data = waal_pca.df,  
  family = gaussian(), 
  chains = 4,          
  cores = 4,           
  iter = 2000          
)

save(rpt_brms.waal, file = "Data_outputs/multivariate_rpt_model_waal.RData")


### Extract variance components ------------------------------------------------

source("RPT_functions.R")

extract_variances.brms(rpt_brms.bba)
extract_variances.brms(rpt_brms.waal)


## Extract BLUPS -------------------------------------------------------------

### Plot the BLUPs alone ----
#BBAL#
blups.bba <- ranef(rpt_brms.bba)$ring
blups_pc1 <-  blups.bba[, "Estimate", "Comp1_Intercept"]
blups_pc1.bba <- data.frame(ring = names(blups_pc1), estimate = as.vector(blups_pc1))
blups_pc2 <-  blups.bba[, "Estimate", "Comp2_Intercept"]
blups_pc2.bba <- data.frame(ring = names(blups_pc2), estimate = as.vector(blups_pc2))

blups_pc1.bba %<>% arrange(estimate) %>% rename(BLUP = estimate)
blups_pc2.bba %<>% arrange(estimate) %>% rename(BLUP = estimate)

blups_pc1.bba$ring <- factor(blups_pc1.bba$ring, levels = blups_pc1.bba$ring) 
blups_pc2.bba$ring <- factor(blups_pc2.bba$ring, levels = blups_pc2.bba$ring) 

pc1_bba.plot <- ggplot(blups_pc1.bba, aes(x = ring, y = BLUP)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BLUPs for PC1", x = "Ring", y = "BLUP") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pc2_bba.plot <- ggplot(blups_pc2.bba, aes(x = ring, y = BLUP)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BLUPs for PC2", x = "Ring", y = "BLUP") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggpubr::ggarrange(pc1_bba.plot, pc2_bba.plot, cols = 2)


#WAAL#
blups.waal <- ranef(rpt_brms.waal)$ring
blups_pc1 <-  blups.waal[, "Estimate", "Comp1_Intercept"]
blups_pc1.waal <- data.frame(ring = names(blups_pc1), estimate = as.vector(blups_pc1))
blups_pc2 <-  blups.waal[, "Estimate", "Comp2_Intercept"]
blups_pc2.waal <- data.frame(ring = names(blups_pc2), estimate = as.vector(blups_pc2))

blups_pc1.waal %<>% arrange(estimate) %>% rename(BLUP = estimate)
blups_pc2.waal %<>% arrange(estimate) %>% rename(BLUP = estimate)

blups_pc1.waal$ring <- factor(blups_pc1.waal$ring, levels = blups_pc1.waal$ring) 
blups_pc2.waal$ring <- factor(blups_pc2.waal$ring, levels = blups_pc2.waal$ring) 

pc1_waal.plot <- ggplot(blups_pc1.waal, aes(x = ring, y = BLUP)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BLUPs for PC1", x = "Ring", y = "BLUP") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pc2_waal.plot <- ggplot(blups_pc2.waal, aes(x = ring, y = BLUP)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BLUPs for PC2", x = "Ring", y = "BLUP") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggpubr::ggarrange(pc1_waal.plot, pc2_waal.plot, cols = 2)


### Plot BLUPS against covariates ----
ind_blups_pc1 <- ranef(rpt_brms.bba)$`ring:season`[, "Estimate", "Comp1_Intercept"]
ind_blups_pc1.bba <- data.frame(ring = names(ind_blups_pc1), estimate = as.vector(ind_blups_pc1))

# Extract the BLUPS
ind_blups_pc1_df.bba <- data.frame(
  ring = sapply(strsplit(ind_blups_pc1.bba$ring, "_"), function(x) x[1]),
  season = sapply(strsplit(ind_blups_pc1.bba$ring, "_"), function(x) x[2]),
  BLUP = ind_blups_pc1.bba$estimate
)

ind_blups_pc1_df.bba %<>% mutate(season = as.numeric(season))

## Merge into metadata
meta.bba <- bbal.df %>% select(c(ring, season, sex, age, colony))
ind_blups_pc1_df.bba <- merge(ind_blups_pc1_df.bba, meta.bba, by = c("ring", "season"), all.x = T) %>%
  distinct()

#### AGE ----

age_plot_pc1.bba <- ggplot(ind_blups_pc1_df.bba, aes(x = age, y = BLUP, col = sex, group = sex)) +
  geom_point() +
  geom_smooth(method = "loess") + 
  scale_colour_manual(values = c("#FFC107", "#085D41"), name = "") +
  labs(x = "Age (years)", y = "BLUP", title = "Black-browed albatross") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95),
        legend.background = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

#### SEX ----

sex_plot_pc1.bba <- ggplot(ind_blups_pc1_df.bba, aes(x = sex, y = BLUP, col = sex)) +
  geom_boxplot() +
  scale_colour_manual(values = c("#FFC107", "#085D41"), name = "") +
  labs(x = "Sex", y = "BLUP", title = "Black-browed albatross") +
  theme_bw() +
  theme(legend.position = "none")


#### INDIVIDUAL ----

my_rings <- unique(ind_blups_pc1_df.bba$ring) 

for (i in 1:length(my_rings)) {
  
  ind_data <- ind_blups_pc1_df.bba %>% 
    filter(ring == my_rings[i]) %>%
    filter(n() > 3)
  
  if(nrow(ind_data) == 0) next
  
  ind_plot <- ind_data %>% ggplot(aes(x = season, y = BLUP)) +
    geom_point() + 
    geom_line(group = 1) +
  theme_bw()
  
  print(ind_plot)

}
  
# Lifetime individual variation
blups_ind <- blups_df.bba %>% group_by(ring) %>%
  summarise(n_recs = n(),
            med_blup = median(BLUP),
            var_blup = sd(BLUP))

# WAAL ----
lmer.waal <- lmer(duration.hrs ~ (1|ring) + (1|season) + (1|ring:season) +
                   phase + sex + age, 
                 data = rep_data %>% filter(species == "waal"))


### Plot the BLUPs alone ----
blups_only.waal <- ranef(lmer.waal)$ring
blups_only.waal$ring <- rownames(blups_only.waal)
blups_only.waal %<>% rename(BLUP = `(Intercept)`) %>% arrange(BLUP)
blups_only.waal$ring <- factor(blups_only.waal$ring, levels = blups_only.waal$ring) 

ggplot(blups_only.waal, aes(x = ring, y = BLUP)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BLUPs for Ring", x = "Ring", y = "BLUP") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

### Plot BLUPS against covariates ----
blups.waal <- ranef(lmer.waal)$`ring:season`

# Extract the BLUPS
blups_df.waal <- data.frame(
  ring = sapply(strsplit(rownames(blups.waal), ":"), function(x) x[1]),
  season = sapply(strsplit(rownames(blups.waal), ":"), function(x) x[2]),
  BLUP = blups.waal$`(Intercept)`
)

blups_df.waal %<>% mutate(season = as.numeric(season))

## Merge into metadata
meta.waal <- rep_data %>% select(c(ring, season, sex, age, colony))
blups_df.waal <- merge(blups_df.waal, meta.waal, by = c("ring", "season"), all.x = T)

#### AGE ----

age_plot.waal <- ggplot(blups_df.waal, aes(x = age, y = BLUP, col = sex, group = sex)) +
  geom_point() +
  geom_smooth(method = "loess") + 
  scale_colour_manual(values = c("#FFC107", "#085D41"), name = "") +
  labs(x = "Age (years)", y = "BLUP", title = "Wandering albatross") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95),
        legend.background = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

#### SEX ----

sex_plot.waal <- ggplot(blups_df.waal, aes(x = sex, y = BLUP, col = sex)) +
  geom_boxplot() +
  scale_colour_manual(values = c("#FFC107", "#085D41"), name = "") +
  labs(x = "Sex", y = "BLUP", title = "Wandering albatross") +
  theme_bw() +
  theme(legend.position = "none")

## COMBINED ----

gridExtra::grid.arrange(age_plot.bba, age_plot.waal)
gridExtra::grid.arrange(sex_plot.bba, sex_plot.waal, ncol = 2)

# ______________________________ ####
# REFERENCE ---------------------------------------------------------------



# BLUPS for Woods Hole team -----------------------------------------------

woods_blups <- rep_waal %>% filter(colony == "crozet" & phase == "incubation" &
                                     !is.na(max_distance.km))

lmer.woods <- lmer(max_distance.km ~ (1|ring) + (1|season) + (1|ring:season) +
                     age, data = rep_waal)


# Extract blups - ring only
blups_only.waal <- ranef(lmer.woods)$ring
blups_only.waal$ring <- rownames(blups_only.waal)
blups_only.waal %<>% rename(BLUP = `(Intercept)`) %>% arrange(BLUP)
blups_only.waal$ring <- factor(blups_only.waal$ring, levels = blups_only.waal$ring) 

write.csv(blups_only.waal, "WAAL_blups.csv", row.names = F)

# Extract the blups - ring + Season
blups.waal <- ranef(lmer.woods)$`ring:season`

# Extract the BLUPS
blups_df.waal <- data.frame(
  ring = sapply(strsplit(rownames(blups.waal), ":"), function(x) x[1]),
  season = sapply(strsplit(rownames(blups.waal), ":"), function(x) x[2]),
  BLUP = blups.waal$`(Intercept)`
)

blups_df.waal %<>% mutate(season = as.numeric(season))

## Merge into metadata
meta.waal <- rep_data %>% select(c(ring, season, age, colony))
blups_df.waal <- merge(blups_df.waal, meta.waal, by = c("ring", "season"), all.x = T)

blups_df.waal %<>% mutate(species = "waal")

write.csv(blups_df.waal, "WAAL_blups.csv", row.names = F)


# LIKELIHOOD RATIOS FOR RANDOM EFFECTS ------------------------------------

random_effects <- c("ring", "ring:season")

# BBAL
LRT_results_list.bba <- lapply(random_effects, function(rndm_eff) {
  lrt_test(formula, rndm_eff, focal_cov, rep_bba.cov, myspecies = "bbal")
})

LRT.bba.df <- do.call(rbind, LRT_results_list.bba)
rownames(LRT.bba.df) <- NULL

# WAAL
LRT_results_list.waal <- lapply(random_effects, function(rndm_eff) {
  lrt_test(formula, rndm_eff, focal_cov, rep_waal.cov, myspecies = "waal")
})

LRT.waal.df <- do.call(rbind, LRT_results_list.waal)
rownames(LRT.waal.df) <- NULL

tik <- Sys.time()
print(paste0("Seconds to complete LRT: ", tik-tok))


rndm_output.full <- do.call("rbind", rndm_output.list)
rndm_output.full %<>% arrange(species) %>%
  mutate(AIC = round(AIC, 2),
         BIC = round(BIC, 2),
         logLik = round(logLik, 2))



### TRIP DURATION - OUT-OF-DATE ----

load("Data_inputs/repeatability_data.RData")

#### Within- vs between- year repeatability ----------------------------------------------------
rpt_tripDur.bbal <- rpt(duration.hrs ~ (1|ring) + (1|season) , 
                        grname = c("ring", "season"),
                        data = rep_data %>% filter(species == "bba"), 
                        datatype = "Gaussian", 
                        nboot = 1000) 

rep_within.bbal <- rpt_tripDur.bbal$R$ring
rep_between.bbal <- rpt_tripDur.bbal$R$season


rep_waal %<>% filter(phase == "incubation" & colony == "crozet" !is.na(total_distance.km)) 

rpt_tripDur.waal <- rpt(total_distance.km ~ (1|ring) + (1|season), 
                        grname = c("ring", "season"),
                        data = rep_waal, 
                        datatype = "Gaussian", 
                        nboot = 1000) 

rep_between.waal <- rpt_tripDur.waal$R$ring



par(mfrow = c(2,2))
plot(rpt_tripDur.waal, grname = "ring", type = "boot", cex.main = 0.8, main = "Within-season repeatability for WAAL")
plot(rpt_tripDur.waal, grname = "season", type = "boot", cex.main = 0.8, main = "Between-season repeatability for WAAL")
plot(rpt_tripDur.bbal, grname = "ring", type = "boot", cex.main = 0.8, main = "Within-season repeatability for BBAL")
plot(rpt_tripDur.bbal, grname = "season", type = "boot", cex.main = 0.8, main = "Between-season repeatability for BBAL")


# Fit models
rpt_tripDur.bbal.cov <- rpt(duration.hrs ~ (1|ring) + (1|season) +
                              phase + sex + age, 
                            grname = c("ring", "season"),
                            data = rep_data %>% filter(species == "bba"), 
                            datatype = "Gaussian", 
                            nboot = 1000) 

rep_within.bbal.cov <- rpt_tripDur.bbal.cov$R$ring
rep_between.bbal.cov <- rpt_tripDur.bbal.cov$R$season

rpt_tripDur.waal.cov <- rpt(duration.hrs ~ (1|ring) + (1|season) +
                              phase + sex + age + colony, 
                            grname = c("ring", "season"),
                            data = rep_data %>% filter(species == "waal"), 
                            datatype = "Gaussian", 
                            nboot = 1000) 

rep_within.waal.cov <- rpt_tripDur.waal.cov$R$ring
rep_between.waal.cov <- rpt_tripDur.waal.cov$R$season


par(mfrow = c(2,2))
plot(rpt_tripDur.waal.cov, grname = "ring", type = "boot", cex.main = 0.8, main = "Within-season repeatability for WAAL")
plot(rpt_tripDur.waal.cov, grname = "season", type = "boot", cex.main = 0.8, main = "Between-season repeatability for WAAL")
plot(rpt_tripDur.bbal.cov, grname = "ring", type = "boot", cex.main = 0.8, main = "Within-season repeatability for BBAL")
plot(rpt_tripDur.bbal.cov, grname = "season", type = "boot", cex.main = 0.8, main = "Between-season repeatability for BBAL")




### Manual repeatabiltiy -------

# Build a basic model

# Starting with a manual LMER approach as I don't quite understand rptR
tripDur_lmm.waal <- lmer(duration.hrs ~ 1 + (1|ring) + (1|season),
                         data = rep_data %>% filter(species == "waal"),
                         REML = T)

## Extract variance components
var_components <- VarCorr(tripDur_lmm.waal)

### Within-year = variance between individuals
within_var <- as.numeric(var_components$ring)
resid_var <- attr(var_components, "sc")^2

rep_within <- within_var / (within_var + resid_var) #  0.04903376

### Between-year = variance within individuals between years
between_var <- as.numeric(var_components$season)

rep_between <- between_var / (between_var + resid_var) #  0.07319174


