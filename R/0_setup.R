
# Create required directories if they don't already exist ----
required_dirs <- c(
  "Output",
  "Output/Figures",
  "Output/Models",
  "Output/Plots",
  "Output/Posteriors",
  "Data",
  "Data/Intermediate",
  "Data/Temporary"
)

invisible(lapply(
  required_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))


# Configuration & options -----

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)

# Ensure explicitly using dplyr select
select <- dplyr::select


# Load custom functions ----

source("R/Functions/RPT_functions.R")
