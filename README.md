# ALB_rpt

Repository for repeatability analyses.

## Scripts

A short description of each script is given below. Briefly, the workflow is as follows:

1.  Process raw GPS/GLS data into labelled trips;
2.  Extract foraging behaviors and create individual trip summaries;
3.  Sample points for habitat analyses and add environmental variables;
4.  Run analyses.

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

## Data inputs

These are the raw and processed datasets used as inputs to the scripts listed above. Files are listed in alphabetical order by their description suffix (i.e. ignoring leading species/colony demarker). Other files are produced during the analyses but are not presented here.

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

-   `ind_meta.RData`: Individual-level metadata, summarised from the demographic dataset (`{BAS/Chize}_all_demo.RData`).

    -   *ring*: Unique bird ID;
    -   *sex*: Sex of the bird;
    -   *species*: Species code;
    -   *birthYr*: Hatch year of the bird.

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
