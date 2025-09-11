# **No individual or heritable signal in foraging of wide-ranging albatrosses**

Natasha Gillies, Richard A. Phillips, Henri Weimerskirch, Jonathan Potts, Denis Réale, Alastair J. Wilson, Frédéric Angelier, Christophe Barbraud, Ashley Bennison, Karine Delord, Prescillia Lemesle, Samuel Peroteau, Andrew G. Wood, José C. Xavier, Samantha C. Patrick.

## Overview

This repository contains scripts and data to recreate the main results and figures of the following in-prep manuscript:

No individual or heritable signal in foraging of wide-ranging albatrosses. Gillies N, Phillips RA, Weimerskirch H, Potts J, Réale D, Wilson AJ, Angelier F, Barbraud C, Bennison A, Delord K, Lemesle P, Peroteau S, Wood AG, Xavier JC, Patrick SC

## Environment

-   R version: 4.5.0 (2025-04-11 ucrt, “How About a Twenty-Six”)

-   Platform: x86_64-w64-mingw32 (Windows)

-   Key R packages required: `adehabitatHR`, `betareg`, `brms`, `corrplot`, `data.table`, `DHARMa`, `dplyr`, `ecmwfr`, `emmeans`, `ggplot2`, `ggpubr`, `ggridges`, `glmmTMB`, `gridExtra`, `igraph`, `kableExtra`, `lubridate`, `MasterBayes`, `magrittr`, `MCMCglmm`, `MASS`, `momentuHMM`, `nadiv`, `pheatmap`, `purrr`, `raster`, `rnaturalearth`, `reticulate`, `sf`, `stringr`, `terra`, `tidyr`, `readxl`

## Scripts

A short description of each script is given below. Briefly, the workflow is as follows:

1.  Process raw GPS/GLS data into labelled trips;
2.  Extract foraging behaviors and create individual trip summaries;
3.  Sample points for habitat analyses and add environmental variables;
4.  Run analyses.

### Functions

-   **RPT_functions.R**: Contains custom functions used throughout the workflow for data processing, behavioural extraction, habitat sampling, and repeatability analyses. This script is sourced by other R and RMarkdown scripts as needed.

### Data processing

These scripts generate the core data files used to run the analyses. They are listed in order of usage.

-   **extract_trips_generalised**: Extracts foraging trips from the original GPS and GLS datasets. This script processes the raw tracking data into labelled trip-level datasets that are used to extract finer-scale behaviour in the 'extract_foraging_behaviours' script.
    -   *Inputs*: `{species}_{colony}_gls_{dates}.RData`; `{species}_{colony}_gls_labelled_{dates}.RData`
    -   *Outputs*: `{species}_{colony}_gps_labelled.RData`; `{species}_{colony}_gls_labelled.RData`
-   **extract_foraging_behaviours**: Uses the labelled GPS and GLS trips to identify behaviours during foraging trips (e.g., wet–dry bouts, speed-based movement bouts). Also merges in demographic and breeding information at the individual and pair level.
    -   *Inputs*: `{species}_{colony}_gps_labelled.RData`; `{species}_{colony}_gls_labelled.RData`; `{BAS/Chize}_demo_complete.RData`
    -   *Outputs*: `{species}_{colony}_individual_trips_summary.RData`
-   **isolate_foraging_sample_points:** Labels GPS tracks with behavioural states, extracts foraging components, and samples available (pseudo-absence) locations for habitat-use analyses.
    -   *Inputs*: `{species}_{colony}_gps_labelled.RData`; `breeding_dates.RData` ; `ind_meta.RData`
    -   *Outputs*: `{species}_{colony}_available_pnts_population_level.RData`
-   **extract_environment**: Extracts environmental variables (sea surface temperature, bathymetry, chlorophyll-a) and matches them to the sampled GPS points for population-level analyses.
    -   *Inputs*: `{species}_{colony}_available_pnts_population_level.RData`
    -   *Outputs*: `{species}_{colony}_sst-bathy-chla_data_population_level.RData`

### Data analyses

These scripts are used to run the main analyses in the text. They are written as RMarkdown files to promote usability and to add important information. They are listed in the order they are presented in the manuscript.

-   **1_RPT_foraging_effort.Rmd**: Analyses for repeatability of foraging effort. Fits models to quantify the consistency of individual foraging behaviour across trips and estimate heritability using the pedigree information.
    -   *Inputs*: `{species}_{colony}_individual_trips_summary.RData`; `ind_meta.RData`; `{species}_{colony}_pedigree.csv`
    -   *Outputs*: Figure 1, Table 1, Figure SX; `behavioural_repeatability_data.RData`; `behavioural_repeatability_results.RData`; `{species}_{colony}_brms_behaviour_repeatability_{behaviour}.RData`; `{species}_{colony}_behaviour_heritability_posteriors.RData`
-   **2_RPT_spatial.Rmd**: Analyses testing for spatial repeatability by examining consistency in foraging bearings and space use, using Bhattacharyya’s affinity.
    -   *Inputs*: `{species}_{colony}_gps_labelled.RData`; `{species}_{colony}_pedigree.csv`
    -   *Outputs*: Figures 2, 3, SX, & Sx; `all_gps.RData`; `{species}_{colony}_ba_comparison.RData`; `ba_means.RData`; `ba_model_summary.RData`; `bearings_data.RData`; `{species}_{colony}_bearing_differences.RData`; `bearing_means.RData`; `bearing_models_summary.RData`
-   **3_RPT_habitat_preference.Rmd**: Analyses for repeatability in habitat use and preference by individual birds. Quantifies individual consistency in habitat selection and links environmental covariates to usage patterns.
    -   *Inputs*: `{species}_{colony}_sst-bathy-chla_data_population_level.RData`; `{species}_{colony}_pedigree.csv`
    -   *Outputs*: Figures 4 & 5, Tables 4 & 5; `{species}_{colony}_glmm_{between/within}_{environment}.RData`; `{species}_{colony}_habitat_data.RData`; `{species}_{colony}_habitat_slopes.RData`; `{species}_{colony}_habitat_used_ranef.RData`; `{species}_{colony}_habitat_used_variances.RData`
-   **4_RPT_heritability.Rmd**: Combines heritability estimates from other analyses, generates summary plots, and compiles a table of results.
    -   *Inputs*: `all_posteriors_spatial_heritability.RData`; `{species}_{colony}_behaviour_heritability_posteriors.RData`; `{species}_{colony}_habitat_pref_heritability_posteriors.RData`
    -   *Outputs*: Figure 6, Table 6

## Data inputs

These are the raw and processed datasets used as inputs to the scripts listed above. Files are listed in alphabetical order by their description suffix (i.e. ignoring leading species/colony demarker). Other files produced during the analyses are not presented here.

-   `{BAS/Chize}_demo_complete.RData`: Demographic and breeding information for individuals and pairs.

    -   *bird1*, *bird2*: Unique IDs for the pair members;
    -   *season*: Breeding season;
    -   *pairID*: Unique identifier for the pair (combination of bird1 and bird2);
    -   *rs*: Reproductive success code;
    -   *colony*: Breeding colony;
    -   *species*: Species code;
    -   *nest*: Nest identifier;
    -   *lay_date*: Date egg was laid;
    -   *hatch_date*: Date chick hatched;
    -   *fail_date*: Date breeding attempt failed (if applicable);
    -   *fledge_date*: Date chick fledged (if successful);
    -   *bird1_sex*, *bird2_sex*: Sex of each bird;
    -   *bird1_birthYr*, *bird2_birthYr*: Hatch year of each bird.

-   `{species}_{colony}_available_pnts_population_level.RData`: Dataset of used and available points for habitat-use analyses, at the population level. Each row represents a used (observed) or available (sampled) location.

    -   *ring*: Unique bird ID;
    -   *season*: Breeding season;
    -   *year*: Calendar year;
    -   *tripID*: Foraging trip identifier;
    -   *date_hourly*: Timestamp rounded to the nearest hour;
    -   *longitude*, *latitude*: Geographic coordinates of the point;
    -   *used*: Indicator (1 = observed point, 0 = sampled available point);
    -   *phase*: Breeding phase during which the point occurred;
    -   *sex*: Sex of the bird;
    -   *dist_col*: Distance from colony (km).

-   `breeding_dates.RData`: Summarised breeding phenology information for individuals, derived from the demographic dataset (`demo_complete`).

    -   *ring*: Unique bird ID;
    -   *season*: Breeding season;
    -   *partner*: Partner bird ID;
    -   *pairID*: Unique pair identifier;
    -   *species*: Species code;
    -   *colony*: Breeding colony;
    -   *lay_date*: Date egg was laid;
    -   *hatch_date*: Date chick hatched.

-   `{species}_{colony}_gls_{dates}.RData`: Raw GLS (geolocator) data, including wet–dry immersion records. Each row represents a recorded timepoint for an individual bird.

    -   *ring*: Unique bird ID;
    -   *season*: Breeding season, defined as the calendar year in which chicks are expected to fledge, i.e. the year following the onset of breeding;
    -   *darvic*: Plastic colour ring code;
    -   *year_started*: First year of tracking;
    -   *datetime*: Timestamp of record;
    -   *immersion*: Immersion value (0-200);
    -   *partner*: Partner bird ID;
    -   *rs*: Reproductive success code.

-   `{species}_{colony}_gps_{dates}.RData`: Raw GPS data from tracking devices. Each row represents a location fix for an individual bird.

    -   *season*: Breeding season during which the track was recorded;
    -   *species*: Species code (BBAL or WAAL);
    -   *ring*: Unique bird ID;
    -   *datetime*: Timestamp of location fix;
    -   *longitude*, *latitude*: Geographic coordinates;
    -   *year*: Calendar year of the fix;
    -   *ringYr*: Bird ID concatenated with year (unique individual-year identifier).

-   `{species}_{colony}_gls_labelled.RData`: GLS trips processed and labelled with trip IDs and behavioural bouts.

    -   *ring*: Unique bird ID;
    -   *season*: Breeding season during which the trip occurred;
    -   *datetime*: Timestamp of record;
    -   *immersion*: Wet–dry immersion value (from GLS logger);
    -   *birdSeason*: Unique identifier combining bird and season;
    -   *dryWet*: Immersion state (dry or wet);
    -   *dur.mins*: Duration of the dry/wet bout in minutes;
    -   *loc*: Location estimate associated with the bout; on foraging trip or at colony (trip or col);
    -   *tripID*: Unique identifier for the foraging trip.

-   `{species}_{colony}_gps_labelled.RData`: GPS trips processed and labelled with trip IDs and behavioural bouts.

    -   *ring*: Unique bird ID;
    -   *season*: Breeding season during which the trip occurred;
    -   *ringYr*: Bird ID concatenated with year (individual–year identifier);
    -   *datetime*: Timestamp of location fix;
    -   *longitude*, *latitude*: Geographic coordinates of fix;
    -   *calc_speed*: Estimated speed (m/s) between consecutive fixes;
    -   *dist_col*: Distance from colony (km);
    -   *loc*: Labelled location category (e.g., colony, sea);
    -   *bout*: Foraging bout identifier;
    -   *boutID*: Unique ID for each bout within a trip;
    -   *duration.mins*, *duration.hours*: Duration of the bout in minutes/hours;
    -   *tripID*: Unique identifier for the foraging trip;
    -   *ID*: Internal unique trip or individual identifier;
    -   *rs*: Reproductive success code for the pair in that season;
    -   *lay_date*: Date egg was laid;
    -   *hatch_date*: Date chick hatched;
    -   *fail_date*: Date breeding attempt failed (if applicable);
    -   *hatch_code*: Hatch outcome code (e.g., success/failure);
    -   *phase*: Breeding phase during which the trip occurred.

-   `{species}_{colony}_habitat_data.RData`: GLS trips processed and labelled with trip IDs and behavioural bouts.

-   `ind_meta.RData`: Individual-level metadata, summarised from the demographic dataset (`{BAS/Chize}_all_demo.RData`).

    -   *ring*: Unique bird ID;
    -   *sex*: Sex of the bird;
    -   *species*: Species code;
    -   *birthYr*: Hatch year of the bird.

-   `{species}_{colony}_pedigree.csv`: Pedigree information for individuals used in repeatability and heritability analyses. When sex was unknown, dams and sires were randomly assigned.

    -   *id*: Unique identifier for each individual bird;
    -   *dam*: Mother’s unique ID;
    -   *sire*: Father’s unique ID;
    -   *Year*: Year of hatching;

-   `{species}_{colony}_sst-bathy-chla_data_population_level.RData`: GPS points with environmental variables appended. Each row corresponds to a used or available location with matched environmental data.

    -   *index*: Internal row identifier;
    -   *ring*: Unique bird ID;
    -   *season*: Breeding season;
    -   *year*: Calendar year;
    -   *tripID*: Foraging trip identifier;
    -   *date_hourly*: Timestamp of the point (rounded to the nearest hour);
    -   *longitude*, *latitude*: Geographic coordinates;
    -   *used*: Indicator for observed (1) vs available (0) point;
    -   *phase*: Breeding phase during which the point occurred;
    -   *sex*: Sex of the bird;
    -   *dist_col*: Distance from colony (km);
    -   *SST*: Sea surface temperature (°C);
    -   *bathy*: Water depth (m);
    -   *rounded_time*: Time of observation rounded for environmental matching;
    -   *chlA*: Chlorophyll-a concentration (mg m⁻³).
