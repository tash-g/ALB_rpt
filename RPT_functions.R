
# * DATA PROCESSING * -----------------------------------------------------

cut_to_breeding <- function(data, medLay, medHatch, hatchTime, broodTime) {
  
  # Process the data
  data %>%
    mutate(rs = ifelse(is.na(rs), "UNKNOWN", rs)) %>%
    mutate(lay_date = ifelse(is.na(lay_date), paste0(season - 1, "-", medLay), as.character(lay_date)),
           hatch_code = ifelse(is.na(hatch_date), "NEW_EST", "MEAS"),
           hatch_date = ifelse(is.na(hatch_date) & rs != "FAILED_EGG", paste0(season, "-", medHatch), as.character(hatch_date)),
           hatch_month = month(as.Date(hatch_date)),
           hatch_date = ifelse(hatch_month > 9 & hatch_code == "NEW_EST", 
                               as.character(as.Date(hatch_date) - 365), as.character(hatch_date))) %>%
    mutate(across(c(lay_date, hatch_date, fail_date), as.Date)) %>%
    select(-c(hatch_month)) %>%
    mutate(flag = ifelse(!is.na(fail_date) & as.Date(datetime) >= fail_date, "flag", "fine")) %>%
    filter(flag == "fine") %>%
    select(-flag) %>%
    mutate(brood_end = ifelse(!is.na(hatch_date), as.character(hatch_date + broodTime),
                              as.character(lay_date + hatchTime))) %>%
    mutate(brood_end = as.Date(brood_end)) %>%
    group_by(ring, season, boutID) %>%
    mutate(start_time = first(datetime, 1), end_time = last(datetime, 1),
           remove = ifelse(as.Date(end_time) < lay_date - 1 | as.Date(end_time) > brood_end, "remove", "keep")) %>%
    filter(remove == "keep") %>%
    select(-remove) %>%
    mutate(phase = ifelse(as.Date(datetime) < hatch_date | is.na(hatch_date), "incubation", "brooding")) %>%
    arrange(ring, datetime) %>%
    select(-brood_end)
}



label_breeding <- function(data, medLay, medHatch) {
  
  # Process the data
  data %>%
    mutate(rs = ifelse(is.na(rs), "UNKNOWN", rs)) %>%
    mutate(lay_date = ifelse(is.na(lay_date), paste0(season - 1, "-", medLay), as.character(lay_date)),
           hatch_code = ifelse(is.na(hatch_date), "NEW_EST", "MEAS"),
           hatch_date = ifelse(is.na(hatch_date) & rs != "FAILED_EGG", paste0(season, "-", medHatch), as.character(hatch_date)),
           hatch_month = month(as.Date(hatch_date)),
           hatch_date = ifelse(hatch_month > 9 & hatch_code == "NEW_EST", 
                               as.character(as.Date(hatch_date) - 365), as.character(hatch_date))) %>%
    mutate(across(c(lay_date, hatch_date, fail_date), as.Date)) %>%
    select(-c(hatch_month)) %>%
    mutate(flag = ifelse(!is.na(fail_date) & as.Date(datetime) >= fail_date, "flag", "fine")) %>%
    filter(flag == "fine") %>%
    select(-flag) %>%
    mutate(phase = ifelse(as.Date(datetime) < hatch_date | is.na(hatch_date), "incubation", "brooding")) %>%
    arrange(ring, datetime) 
}


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



# * HERITABILITY FUNCTIONS * ----------------------------------------------

get_ancestors <- function(ids, pedigree) {
  all_ids <- ids
  new_ids <- ids
  
  repeat {
    # Get dams and sires of current new_ids
    parents <- pedigree[pedigree$id %in% new_ids, c("dam", "sire")]
    parents <- unique(na.omit(unlist(parents)))
    
    # Stop if no new parents found
    parents <- setdiff(parents, all_ids)
    if (length(parents) == 0) break
    
    # Add new parents to the list
    all_ids <- c(all_ids, parents)
    new_ids <- parents
  }
  
  return(unique(all_ids))
}


# * MODEL DIAGNOSES * -----------------------------------------------------

get_auc <- function(model) {
  model_data <- model.frame(model)
  response <- model_data$used
  predicted_probs <- predict(model, type = "response")
  roc_obj <- roc(response, predicted_probs)
  auc_value <- auc(roc_obj)
  return(cat("AUC:", auc_value, "\n"))
}



overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  chisq <- sum(rp^2)
  c(chisq = chisq, ratio = chisq / rdf, rdf = rdf)
}


# * PLOTTING FUNCTIONS * --------------------------------------------------

make_plot <- function(df) {
  df <- df %>%
    mutate(ring_ordered = forcats::fct_reorder(ID, estimate))
  
  ggplot(df, aes(x = estimate, y = ring_ordered, col = env)) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    scale_color_manual(values = env_colours, guide = "none") +
    facet_wrap(~group, scales = "free") +
    labs(title = unique(df$env),
         x = "Random intercept estimate",
         y = "ID") +
    theme_bw() +
    theme(axis.text.y = element_blank())}


# * REPEATABILITY CALCULATIONS * ------------------------------------------

calc_rpt_brms <- function(brm_model, data, x_mean) {
  # Extract posterior draws as dataframe
  posterior <- as_draws_df(brm_model)
  
  # Extract SDs and correlations from posterior
  sd_ring_season <- posterior$sd_ring_season__Intercept
  sd_slope <- posterior$sd_ring_season__time_since_first
  cor_is <- posterior$cor_ring_season__Intercept__time_since_first
  
  var_ring_season_intercept <- sd_ring_season^2
  var_slope <- sd_slope^2
  var_ring <- posterior$sd_ring__Intercept^2
  var_season <- posterior$sd_season__Intercept^2
  var_resid <- posterior$sigma^2
  
  # Calculate ring_season variance including slope variance & covariance
  var_ringSeason_x <- var_ring_season_intercept +
    (x_mean^2) * var_slope +
    2 * x_mean * cor_is * sd_ring_season * sd_slope
  
  ## Within-season variance includes ID variance
  var_ID_within_season <- var_ring + var_ringSeason_x
  
  # Total variance
  total_var <- var_ringSeason_x + var_ring + var_season + var_resid
  
  # Repeatabilities
  rpt_within <- var_ID_within_season / total_var
  rpt_across <- var_ring / total_var
  
  # Combine into tibble/data.frame with all posterior draws
  result <- tibble(
    var_ringSeason_x = var_ID_within_season,
    var_ring = var_ring,
    var_season = var_season,
    var_resid = var_resid,
    sd_slope = sd_slope,
    rpt_within = rpt_within,
    rpt_across = rpt_across
  )
  
  return(result)
}




get_rnd_effects <- function(mod_fit, hab_var) {
  
  random_effects <- ranef(mod_fit, condVar = TRUE)[[1]][1]
  slopes <- random_effects[[names(random_effects)[1]]][[hab_var]]
  intercepts <- random_effects[[names(random_effects)[1]]][["(Intercept)"]]        
  
  rnd_eff.df <- data.frame(cbind(intercepts, slopes))
  colnames(rnd_eff.df) <- c("intercept", "slope")
  rnd_eff.df$ID = rownames(random_effects[[1]])
  rownames(rnd_eff.df) <- NULL
  
  return(rnd_eff.df)
}



# * SPATIAL FUNCTIONS * ---------------------------------------------------

angle_diff <- function(b1, b2) {
  abs_diff <- abs(b1 - b2) # %% 360
  return(ifelse(abs_diff > 180, 360 - abs_diff, abs_diff))
}



convert_to_spdf <- function(trip, projection) {
  
  mytrip <- trip %>% ungroup() %>% select(c(longitude, latitude))
  coordinates(mytrip) <- ~longitude+latitude
  proj4string(mytrip) <- proj.dec
  return(mytrip)
}




calculate_bearing_diffs <- function(df, bearing_col_name, comparison_type = c("within", "between"), output_col_name, trip_starts) {
  
  comparison_type <- match.arg(comparison_type)
  bearing_col_sym <- rlang::sym(bearing_col_name)
  
  # Select relevant columns
  df <- df %>%
    ungroup() %>%
    select(ring, boutID, season, phase, !!bearing_col_sym)
  
  if (comparison_type == "within") {
    
    out <- df %>%
      group_by(ring) %>%
      filter(n() >= 2) %>%
      summarise(
        trip_comparisons = list(combn(boutID, 2, paste, collapse = " vs ")),
        bearing_diff_tmp = list(combn(!!bearing_col_sym, 2, function(x) angle_diff(x[1], x[2]))),
        .groups = "drop" ) %>%
      unnest(cols = c(trip_comparisons, bearing_diff_tmp)) %>%
      rename(!!output_col_name := bearing_diff_tmp) %>%
      separate(trip_comparisons, into = c("tripID_1", "tripID_2"), sep = " vs ") %>%
      left_join(df, by = c("ring", "tripID_1" = "boutID")) %>%
      rename(season_1 = season, phase_1 = phase) %>%
      left_join(df, by = c("ring", "tripID_2" = "boutID")) %>%
      rename(season_2 = season, phase_2 = phase) %>%
      mutate(
        comp_ind = "within_ring",
        comp_season = ifelse(season_1 == season_2, "within_season", "between_season"),
        comp_phase = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")
      ) %>%
      filter(comp_phase == "same_phase") %>%
      mutate(ring_2 = ring) %>%
      rename(ring_1 = ring) %>%
      select(ring_1, tripID_1, season_1, ring_2, tripID_2, season_2, 
             !!output_col_name, comp_ind, comp_season)
    
    out %<>%
      left_join(trip_starts, by = c("ring_1" = "ring", "tripID_1" = "tripID")) %>%
      rename(start_time_1 = start_time) %>%
      left_join(trip_starts, by = c("ring_2" = "ring", "tripID_2" = "tripID")) %>%
      rename(start_time_2 = start_time) %>%
      mutate(inter_obs_dist = abs(as.numeric(difftime(start_time_1, start_time_2, units = "days"))))
    
  } else if (comparison_type == "between") {
    
    df1 <- df %>%
      rename_with(~ paste0(.x, "_1"))
    
    df2 <- df %>%
      rename_with(~ paste0(.x, "_2"))
    
    out <- merge(df1, df2, by = NULL) %>%
      filter(ring_1 != ring_2, phase_1 == phase_2) %>%
      mutate(
        !!output_col_name := angle_diff(!!rlang::sym(paste0(bearing_col_name, "_1")),
                                        !!rlang::sym(paste0(bearing_col_name, "_2"))),
        comp_ind = "between_ring",
        comp_season = ifelse(season_1 == season_2, "within_season", "between_season") ) %>%
      mutate(
        tripID_1.tmp = paste0(ring_1, boutID_1),
        tripID_2.tmp = paste0(ring_2, boutID_2),
        comp_ID = paste0(pmin(tripID_1.tmp, tripID_2.tmp), "_", pmax(tripID_1.tmp, tripID_2.tmp)) ) %>%
      distinct(comp_ID, .keep_all = TRUE) %>%
      select(ring_1, boutID_1, season_1,
             ring_2, boutID_2, season_2,
             !!output_col_name, comp_ind, comp_season) %>%
      rename(tripID_1 = boutID_1, tripID_2 = boutID_2)
    
    out %<>%
      left_join(trip_starts, by = c("ring_1" = "ring", "tripID_1" = "tripID")) %>%
      rename(start_time_1 = start_time) %>%
      left_join(trip_starts, by = c("ring_2" = "ring", "tripID_2" = "tripID")) %>%
      rename(start_time_2 = start_time) %>%
      mutate(inter_obs_dist = abs(as.numeric(difftime(start_time_1, start_time_2, units = "days"))))
    
  }
  
  return(out)
}





calc_speed <- function(mydf) {
  
  n = nrow(mydf)
  
  mydf$dist_next <- c(gcd.hf(mydf$longitude[2:n],
                             mydf$latitude[2:n],
                             mydf$longitude[1:(n-1)],
                             mydf$latitude[1:(n-1)]),NA)
  
  dt = c(as.numeric(difftime(mydf$datetime[2:n], mydf$datetime[1:(n-1)], 
                             units = "secs")), NA)
  mydf$calc_speed <- mydf$dist_next*1000/dt
  
  return(mydf$calc_speed)
  
}




gcd.hf <- function(long1, lat1, long2, lat2) { 
  R <- 6371 # Earth mean radius [km]
  deg2rad <- function(deg) return(deg*pi/180)
  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  
  coo2 <- function(x){2 * asin(min(1,sqrt(x)))}
  c <- lapply(a, coo2)
  c <- do.call("c", c)
  d = R * c
  return(d) # Distance in km
}




