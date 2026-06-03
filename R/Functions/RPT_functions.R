
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



format_numbers <- function(x) {
  ifelse(
    abs(x) < 0.0001,
    "<0.0001",
    ifelse(
      abs(x) >= 0.01,
      format(round(x, 2), nsmall = 2),  # two decimal places
      signif(x, 2)                       # 2 significant figures
    )
  )
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


scale_2SD <- function(variable) { (variable - mean(variable, na.rm = TRUE)) / sd(variable, na.rm = TRUE) * 2 }


# * HABITAT PREFERENCE FUNCTIONS * ----------------------------------------

# model = mod_within
# env_name = env_vars[i]
# levelX = "within"
# n_sim = n_sim
# envX = env_current
# df_val = df_ID_season


bootstrap_condR <- function(model, env_name, levelX, n_sim, envX, speciesX) {
  
  habitat_var <- paste0(env_name, "_scaled")
  if(speciesX == "WAAL" & env_name == "chl") { habitat_var = "chl_scaled_log"}
  xrange_90 <- quantile(envX, probs = c(0.05, 0.95), na.rm = TRUE)
  
  if(levelX == "between") {
    grp_name <- "ID"
    df_val   <- df_ID
  } else {
    grp_name <- "ID_season"
    df_val   <- df_ID_season  }
  
  # Extract fixed effects and their covariance matrix
  beta_hat <- fixef(model)$cond
  vcov_beta <- vcov(model)$cond
  vc <- VarCorr(model)$cond
  
  # Extract variance components (point estimates)
  var_int <- as.numeric(vc[[grp_name]][1,1])
  grp_slope_name <- paste0(grp_name, ".1")
  var_slope <- as.numeric(vc[[grp_slope_name]][1,1])
  var_tripID <- as.numeric(vc$tripID[1,1])
  var_resid <- (pi^2)/3
  
  # Compute condR 
  condR_sim <- list()
  
  for (n in 1:n_sim ) {
    
    print(paste0("Iteration ", n, " of ", n_sim))
    
    ## 1. Sample fixed effects from multivariate normal
    beta_sim <- mvrnorm(1, mu = beta_hat, Sigma = vcov_beta)
    
    ## 2. Sample variance components from scaled chi-square approximation
    var_int_sim <- var_int * rchisq(1, df = df_val) / df_val
    var_slope_sim <- var_slope * rchisq(1, df = df_val) / df_val
    var_tripID_sim <- var_tripID * rchisq(1, df = df_tripID) / df_tripID
    
    ## 3. Make the covariance matrix
    Sigma_sim = matrix(c(var_int_sim, 0, 0, var_slope_sim), nrow = 2)
    
    ## 4. Compute conditional repeatability
    sim_out <- condR(
      betaX = beta_sim[habitat_var], 
      meanX = mean(envX), 
      Vx = var(envX), 
      Vr = var_resid, 
      Sigma = Sigma_sim, 
      Vo = var_tripID_sim, 
      xrange = xrange_90)
    
    condR_sim[[n]] <- sim_out$condR %>% mutate(sim = n)
    
  }
  
  # Summarise the output
  bind_rows(condR_sim) %>%
    group_by(x) %>%
    summarise(median = median(condR, na.rm = TRUE),
              lower = quantile(condR, 0.025, na.rm = TRUE),
              upper = quantile(condR, 0.975, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(level = levelX)
}


# env_vars = env_vars.pref
# spec_col = "bbal_birdis"
# meta_df = pnts_meta

extract_habitat_slopes <- function(env_vars, spec_col, meta_df) {
  
  # Map over the environmental variables to create a list of dataframes
  slopes_list <- lapply(env_vars, function(env) {
    
    # Load model
    file_path <- paste0("Data_outputs/", spec_col, "_glmm_within_", env, ".RData")
    load(file_path)
    
    # Extract the model object 
    obj_name <- paste0(spec_col, "_glmm.within_", env)
    my_mod <- get(obj_name)
    
    # Ger random slopes
    slopes <- ranef(my_mod)$cond$tripID
    
    # Make the dataframe
    df <- data.frame(
      tripID = rownames(slopes),
      slope = slopes[, 2],
      env = env) %>%
      left_join(meta_df, by = "tripID") %>%
      separate(tripID, into = c("ring", "season", "tripID"), sep = "_", remove = FALSE) %>%
      relocate(ring, season, tripID, phase, env, slope)
    
    return(df)
  })
  
  # Combine everything
  final_df <- bind_rows(slopes_list)
  return(final_df)
}



# model = mod_within
# env_name = env_vars[i]
# levelX = "within"
# envX = env_current
# speciesX = "WAAL"
# colonyX = "Bird Island"

extract_pref_metrics <- function(model, env_name, levelX, envX, speciesX, colonyX) {
  
  habitat_var <- paste0(env_name, "_scaled")
  if(speciesX == "WAAL" & env_name == "chl") { habitat_var = "chl_scaled_log"}
  
  grp_name = ifelse(levelX == "between", "ID", "ID_season")
  xrange_90 <- quantile(envX, probs = c(0.05, 0.95), na.rm = TRUE)
  
  # Extract slope variance + CIs 
  mod_CI <- confint(model, component = "cond", method = "Wald")
  target_row <- grep(paste0("Std.Dev.*\\|", grp_name, ".1"), rownames(mod_CI), value = TRUE)
  
  s_var <- mod_CI[target_row, "Estimate"]^2
  s_LCL <- mod_CI[target_row, "2.5 %"]^2
  s_UCL <- mod_CI[target_row, "97.5 %"]^2
  
  # Estimate conditional repeatability
  vc <- VarCorr(model)$cond
  my_sigma <- matrix(c(as.numeric(vc[[grp_name]]), 0, 0, as.numeric(vc[[paste0(grp_name, ".1")]])), nrow = 2)
  
  res <- condR(betaX = fixef(model)$cond[habitat_var], meanX = mean(envX), Vx = var(envX), 
               Vr = (pi^2)/3, Sigma = my_sigma, Vo = as.numeric(vc$tripID), xrange = xrange_90)
  
  df <- data.frame(species = speciesX, colony = colonyX, level = levelX, env = env_name,
             slope_var = s_var, LCL_slope_var = s_LCL, UCL_slope_var = s_UCL,
             R2 = res$R2s, marginal_R2 = res$Rmar,
             med_rpt = median(res$condR$condR),
             LCL_rpt = quantile(res$condR$condR, 0.025, na.rm = TRUE),
             UCL_rpt = quantile(res$condR$condR, 0.975, na.rm = TRUE)) %>%
    relocate(c(species, colony, env, level), .before = slope_var)
  
  rownames(df) <- NULL
  
  return(df)
}



#df <- split_df[[1]]

make_plot <- function(df) {
  
  df %<>% group_by(env, group) %>%
    arrange(estimate) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  ggplot(df, aes(x = estimate, y = rank, col = env)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    scale_color_manual(values = env_colours, guide = "none") +
    facet_wrap(~group, scales = "free") +
    labs(title = unique(df$env),
         x = "Random intercept estimate",
         y = "ID") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

}



# data = waal_birdis_pnts_all
# env_vars = env_vars.used
# speciesX = "WAAL"
# colonyX = "Bird Island"
# nsim = 20
# mod_form = "~ dist_col + phase + (1| season) + (1 | tripID) + (1 | ID_season) "

model_habitat_use <- function(data, env_vars, speciesX, colonyX, nsim = 100, mod_form) {
  
  var_results <- list()
  ranef_results <- list()
  
  for (env in env_vars) {
    
    message(paste0("Modelling ", env, " for ", speciesX, " ", colonyX, "..."))
    
    # Fit Model
    form <- as.formula(paste(env, mod_form))
    mod <- lmer(form, data = data)
    
    # Extract variance components & calc repeatability
    vc <- as.data.frame(VarCorr(mod))
    v_total  <- sum(vc$vcov)
    
    # Bootstrap CIs
    set.seed(123)
    
    ## Adjust for WAAL Bird Is (no across-season)
    
    if (speciesX == "WAAL" & colonyX == "Bird Island") {
      
      boot <- bootMer(mod, FUN = repeat_func_waBI, nsim = nsim, .progress = "txt")
      rpt_CI <- apply(boot$t, 2, function(x) quantile(x, c(0.025, 0.975)))
      
      v_within <- vc[vc$grp == "ID_season", "vcov"]
      v_across <- NA
      
      rpt_within <- v_within / v_total
      rpt_across <- NA
      
      LCL_w <- rpt_CI[1, 1]; UCL_w <- rpt_CI[2, 1]
      LCL_a <- NA; UCL_a <- NA
      
    } else {
      
      boot <- bootMer(mod, FUN = repeat_func, nsim = nsim, .progress = "txt")
      rpt_CI <- apply(boot$t, 2, function(x) quantile(x, c(0.025, 0.975)))
      
      v_across <- vc[vc$grp == "ID", "vcov"]
      v_within <- vc[vc$grp == "ID_season", "vcov"] + v_across
      
      rpt_across <- v_across / v_total
      rpt_within <- v_within / v_total
      
      LCL_w <- rpt_CI[1, 1]; UCL_w <- rpt_CI[2, 1]
      LCL_a <- rpt_CI[1, 2]; UCL_a <- rpt_CI[2, 2]
      
    }
    
  
    # Store metrics
    var_results[[env]] <- vc %>%
      rename(Effect = grp, Variance = vcov, `Std. dev` = sdcor) %>%
      mutate(
        Environment = env,
        Species = speciesX,
        Colony = colonyX,
        rpt_within = v_within / v_total,
        rpt_across = v_across / v_total,
        LCL_within = LCL_w, UCL_within = UCL_w,
        LCL_across = LCL_a, UCL_across = UCL_a) %>%
      relocate(Species, Colony, Environment, Effect, Variance, `Std. dev`, rpt_across, LCL_across,
               UCL_across, rpt_within, LCL_within, UCL_within) %>%
      select(-c(var1, var2))
    
    # Store random effects
    ranef_results[[env]] <- broom.mixed::tidy(mod, effects = "ran_vals", conf.int = TRUE) %>%
      #filter(group == "tripID") %>%
      mutate(env = env, species = speciesX, colony = colonyX) %>%
      #separate(level, into = c("ring", "season", "trip")) %>%
      relocate(c(species, colony), .before = effect)  
    }
  
  return(list(
    variance = bind_rows(var_results),
    ranefs = bind_rows(ranef_results)
  ))
}



prep_used_data <- function(data, thin_hours) {
  data %>%
    filter(used == 1) %>%
    mutate(time_thin = floor_date(date_hourly, thin_hours)) %>%
    group_by(tripID, time_thin) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(tripID, date_hourly)
}


process_categorical_slopes <- function(model, levelX, n_sim) {
  
  if(levelX == "between") {
    grp_name <- "ID"
    df_val   <- df_ID
  } else {
    grp_name <- "ID_season"
    df_val   <- df_ID_season  }
  
  # Extract fixed effects and their covariance matrix
  vc <- VarCorr(model)$cond
  V_int <- as.numeric(vc[[grp_name]][1,1])
  V_tripID <- as.numeric(vc$tripID[1,1])
  
  # Extract variance components (point estimates)
  grp_slope_name <- paste0(grp_name, ".1")
  V_slopes_matrix <- vc[[grp_slope_name]]
  V_slope_vec <- diag(V_slopes_matrix)
  V_resid <- (pi^2) / 3
  
  # Get CIs for slope variance
  mod_CI <- confint(model, component = "cond", method = "Wald")
  target_rows <- grep(paste0("Std.Dev.*\\|", grp_slope_name), rownames(mod_CI), value = TRUE)
  
  # Extract and square to get variance CIs
  slope_CI_df <- data.frame(
    env = gsub("Std.Dev.|[|].*", "", target_rows), # Clean names
    slope_var = mod_CI[target_rows, "Estimate"]^2,
    LCL_slope_var = mod_CI[target_rows, "2.5 %"]^2,
    UCL_slope_var = mod_CI[target_rows, "97.5 %"]^2  )
  
  # Get repeatability uncertainty
  rpt_sim <- matrix(NA, nrow = n_sim, ncol = length(V_slope_vec))
  
  for(i in seq_len(n_sim)) {
    var_int_sim <- V_int * rchisq(1, df_val) / df_val
    var_tripID_sim <- V_tripID * rchisq(1, df_tripID) / df_tripID
    
    for(j in seq_along(V_slope_vec)) {
      var_slope_sim <- V_slope_vec[j] * rchisq(1, df_val) / df_val
      V_rand_sim <- var_int_sim + var_slope_sim + var_tripID_sim
      rpt_sim[i, j] <- V_rand_sim / (V_rand_sim + V_resid)
    }
  }
  
  # Summarise
  rpt_summary <- data.frame(
    env = names(V_slope_vec),
    R2 = NA,
    marginal_R2 = NA,
    med_rpt = apply(rpt_sim, 2, median),
    LCL_rpt = apply(rpt_sim, 2, quantile, probs = 0.025),
    UCL_rpt = apply(rpt_sim, 2, quantile, probs = 0.975))
  
  final_df <- slope_CI_df %>%
    left_join(rpt_summary, by = "env") %>%
    mutate(species = "BBAL", colony = "Kerguelen", level = levelX) %>%
    relocate(species, colony, env, level)
  
  return(final_df)
}



process_habitat_data <- function(df, avail_ratio, min_recs_per_trip) {
  
  require(lubridate)
  
  df_clean <- df %>%
    # Remove NAs and set distance to km; remove birds with too few data points (i.e. fewer than a day's worth of data)
    filter(!is.na(SST), !is.na(used), !is.na(bathy), !is.na(chlA), !is.na(phase)) %>%
    mutate(dist_col.km = dist_col / 1000,
           ID = ring,
           season = ifelse(as.numeric(month(date_hourly)) < 9, 
                      as.numeric(year(date_hourly)), 
                      as.numeric(year(date_hourly)) + 1)) %>%
    mutate(raw_trip = gsub("_", "", tripID),
           ID_season = paste(ID, season, sep = "_"),
           tripID = paste(ID, season, raw_trip, sep = "_")) %>%
    select(-c(index, raw_trip, sex)) %>%
    
    # Resample to maintain used:available ratio
    group_by(tripID, date_hourly) %>%
    group_modify(~ {
      used_pts <- .x %>% filter(used == 1)
      avail_pts <- .x %>% filter(used == 0)
      target_avail <- nrow(used_pts) * avail_ratio
      
      if (nrow(avail_pts) >= target_avail) {
        avail_pts <- avail_pts %>% slice_sample(n = target_avail)
        bind_rows(used_pts, avail_pts)
      } else {
        tibble() # Drop hours with insufficient background points
      }
    }) %>%
    ungroup() %>%
    
    # Drop individuals with insufficient data
    group_by(tripID) %>%
    filter(n() >= min_recs_per_trip) %>% # Minimum points for a valid trip
    ungroup() %>%
    
    group_by(ID) %>%
    filter(n_distinct(tripID) >= 2) %>% # MUST have at least 2 trips for repeatability
    ungroup() %>%
    
    # Scale environmental variables 
    mutate(
      sst_scaled      = scale_2SD(SST),
      bathy_scaled    = scale_2SD(bathy),
      chl_scaled      = scale_2SD(chlA),
      dist_col_scaled = scale_2SD(dist_col) )
  
  return(df_clean)
}

# * HERITABILITY FUNCTIONS * ----------------------------------------------

extractVar <- function(model_draws, param_name) {
  
  V_cos = model_draws[[paste0("sd_", param_name, "__cosb_Intercept")]]^2
  V_sin = model_draws[[paste0("sd_", param_name, "__sinb_Intercept")]]^2
  
  V_out = V_cos + V_sin
}


# type = "slope"
# env_name = "bathy"
# data = bbal_birdis_pref_h2.df
# Amat = bbal_birdis_pref_h2.Amat
# priors_animal = priors_animal
# output_path = "Data_outputs/bbal_birdis_"

fit_habitat_heritability_model <- function(type, env_name, data, Amat, priors, output_path) {
  
  # Construct the formula
  formula_str <- paste0(env_name, "_", type," ~ 1 + (1 | ring_pe) + (1 | gr(animal, cov = Amat))")
  model_formula <- as.formula(formula_str)
  
  message("Fitting model for: ", env_name, "...")
  
  model <- brm(
    formula = model_formula,
    prior = priors,
    data = data,
    data2 = list(Amat = Amat),
    chains = 4, cores = 4, iter = 8000, warmup = 4000,
    control = list(max_treedepth = 15, adapt_delta = 0.99))
  
  # Save using the environment name in the filename
  label = ifelse(type == "estimate", "used", "pref")
  saveRDS(model, file = paste0(output_path, "h2_", label, "_", tolower(env_name), ".rds"))
  
  return(model)
}




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


get_h2 <- function(model) {
  
  v_animal <- (VarCorr(model, summary = FALSE)$animal$sd)^2
  v_pe <- (VarCorr(model, summary = FALSE)$ring_pe$sd)^2
  v_r <- (VarCorr(model, summary = FALSE)$residual$sd)^2
  
  h2_est <- v_animal[, 1] / (v_animal[, 1] + v_pe[, 1] + v_r[, 1])
  
  return(h2_est)
  
}


prep_used_heritability_data <- function(ranef_data, pedigree_data) {
  
  # Get one row per trip, with columns for each env estimate
  wide_data <- ranef_data %>%
    filter(group == "tripID") %>%
    separate(level, into = c("ring", "season", "trip"), extra = "drop") %>%
    select(ring, season, trip, estimate, env) %>%
    pivot_wider(
      names_from = env, 
      values_from = estimate,
      names_glue = "{env}_estimate") %>%
    mutate(ring_pe = ring, animal = ring)
  
  # Clean Pedigree
  ped <- pedigree_data %>%
    mutate(across(c(id, dam, sire), as.character))
  
  # Create Amat
  Amat <- make_Amat(wide_data, ped)
  
  return(list(data = wide_data, Amat = Amat))
}



prep_slope_heritability_data <- function(slope_data, pedigree_data) {
  
  wide_data <- slope_data %>%
    pivot_wider(
      names_from = env, 
      values_from = slope,
      names_glue = "{env}_slope") %>%
    mutate(ring_pe = ring, animal = ring)
  
  # Filter for Repeatability
  wide_data <- wide_data %>% 
    group_by(ring) %>% 
    filter(n_distinct(tripID) > 1) %>% 
    ungroup()
  
  # Clean Pedigree
  ped <- pedigree_data %>%
    mutate(across(c(id, dam, sire), as.character))
  
  # Create Amat
  Amat <- make_Amat(wide_data, ped)
  
  return(list(data = wide_data, Amat = Amat))
}



prune_pedigree <- function(focal_data, pedigree_df) {
  
  focal_data$ring <- as.character(focal_data$ring)
  
  # Extract IDs from model data and pedigree
  model_ids <- unique(focal_data$ring)
  
  # First prune pedigree
  ped_ids_needed <- get_ancestors(model_ids, pedigree_df)
  pedigree_pruned <- pedigree_df[pedigree_df$id %in% ped_ids_needed, ]
  
  # Add missing IDs as founders
  missing_ids <- setdiff(model_ids, pedigree_pruned$id)
  if(length(missing_ids) > 0) {
    founders <- data.frame(id = missing_ids, dam = NA, sire = NA)
    pedigree_extended <- rbind(pedigree_df, founders)
  } else {
    pedigree_extended <- pedigree_df
  }
  
  # Prune again with extended pedigree
  ped_ids_needed <- get_ancestors(model_ids, pedigree_extended)
  pedigree_pruned <- pedigree_extended[pedigree_extended$id %in% ped_ids_needed, ]
  
  return(pedigree_pruned)
}


# focal_data = mod_data
# pedigree_df = pedigree.df

make_Amat <- function(focal_data, pedigree_df) {

  focal_data$ring <- as.character(focal_data$ring)
  model_ids <- unique(focal_data$ring)
  
  ped <- pedigree_df[,c("id", "dam", "sire")]
  
  # Add missing IDs as founders
  missing_ids <- setdiff(model_ids, ped$id)
  
  if(length(missing_ids) > 0) {
    founders <- data.frame(id = missing_ids, dam = NA, sire = NA)
    ped <- rbind(ped, founders)
  } 
  
  # Prune pedigree
  ped_ids_needed <- get_ancestors(model_ids, ped)
  ped <- ped[ped$id %in% ped_ids_needed, ]
  
  #intersect(unique(ped$dam), unique(ped$sire))
  
  # Build the additive genetic relationship matrix
  names(ped) <- c("animal", "dam", "sire")
  ped <- MasterBayes::orderPed(ped)
  Amat <- as.matrix(nadiv::makeA(ped))

  return(Amat)
}




simulate_pheno <- function(ped, VA, VE) {
  # Get relationship matrix A from the pedigree
  A <- as.matrix(nadiv::makeA(ped))
  A <- A / mean(diag(A))
  
  # Simulate additive genetic values
  a <- mvrnorm(n = 1, mu = rep(0, nrow(A)), Sigma = VA * A)
  
  # Simulate environmental noise
  e <- rnorm(nrow(A), mean = 0, sd = sqrt(VE))
  
  # Phenotype = genetic + environmental
  pheno <- a + e
  
  data.frame(animal = ped$animal, trait_1 = pheno)
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



# * REPEATABILITY CALCULATIONS * ------------------------------------------

calc_rpt_brms_VM <- function(brm_model) {
  # Extract posterior draws
  posterior <- as_draws_df(brm_model)
  
  # Extract SDs of random effects and kappa
  sd_ring         <- posterior$sd_ring__Intercept
  sd_ring_season  <- posterior$sd_ring_season__Intercept
  sd_season       <- posterior$sd_season__Intercept
  kappa           <- posterior$kappa
  
  # Approximate variance on circular scale
  var_ring        <- 2 * (1 - exp(-0.5 * sd_ring^2))
  var_ring_season <- 2 * (1 - exp(-0.5 * sd_ring_season^2))
  var_season      <- 2 * (1 - exp(-0.5 * sd_season^2))
  var_resid       <- 2 * (1 - exp(-0.5 / kappa))
  
  # Within-season variance (ring + ring_season)
  var_within <- var_ring + var_ring_season
  
  # Total variance
  total_var <- var_ring + var_ring_season + var_season + var_resid
  
  # Repeatabilities
  rpt_within  <- var_within / total_var
  rpt_across  <- var_ring / total_var
  
  # Combine into tibble with posterior draws
  result <- tibble(
    var_ring        = var_ring,
    var_ring_season = var_ring_season,
    var_season      = var_season,
    var_resid       = var_resid,
    rpt_within      = rpt_within,
    rpt_across      = rpt_across
  )
  
  return(result)
}


calc_rpt_brms_VM_waBI <- function(brm_model) {
  # Extract posterior draws
  posterior <- as_draws_df(brm_model)
  
  # Extract SDs of random effects and kappa
  sd_ring_season  <- posterior$sd_ring_season__Intercept
  sd_season       <- posterior$sd_season__Intercept
  kappa           <- posterior$kappa
  
  # Approximate variance on circular scale
  var_ring_season <- 2 * (1 - exp(-0.5 * sd_ring_season^2))
  var_season      <- 2 * (1 - exp(-0.5 * sd_season^2))
  var_resid       <- 2 * (1 - exp(-0.5 / kappa))
  
  # Within-season variance (ring + ring_season)
  var_within <- var_ring_season
  
  # Total variance
  total_var <- var_ring_season + var_season + var_resid
  
  # Repeatabilities
  rpt_within  <- var_within / total_var

  # Combine into tibble with posterior draws
  result <- tibble(
    var_ring        = NA,
    var_ring_season = var_ring_season,
    var_season      = var_season,
    var_resid       = var_resid,
    rpt_within      = rpt_within,
    rpt_across      = NA
  )
  
  return(result)
}




calc_rpt_brms_slope <- function(brm_model, data, x_mean) {
  # Extract posterior draws as dataframe
  posterior <- as_draws_df(brm_model)
  
  # Variance components
  var_ring <- posterior$sd_ring__Intercept^2
  var_ringSeason <- posterior$sd_ring_season__Intercept^2
  var_season <- posterior$sd_season__Intercept^2
  var_resid <- posterior$sigma^2
  
  # Individual variance within season
  var_ID_within_season <- var_ring + var_ringSeason
  
  # Total variance
  total_var <- var_ring + var_ringSeason + var_season + var_resid
  
  # Repeatabilities
  rpt_within <- var_ID_within_season / total_var
  rpt_across <- var_ring / total_var
  
  result <- tibble(
    var_ringSeason = var_ringSeason,
    var_ring = var_ring,
    var_season = var_season,
    var_resid = var_resid,
    rpt_within = rpt_within,
    rpt_across = rpt_across
  )
  
  return(result)
}


extract_output.lm <- function(model) {
  
  s <- summary(model)
  cf <- s$coefficients$mean
  ci <- confint(model, level = 0.95)[1:nrow(cf),]
  
  data.frame(
    Term      = rownames(cf),
    Estimate  = cf[, 1],
    Std_Error = cf[, 2],
    t_value   = cf[, 3],
    p_value   = cf[, 4],
    CI_lower  = ci[, 1],
    CI_upper  = ci[, 2],
    row.names = NULL
  )
}


get_rnd_effects <- function(mod_fit, hab_var) {
  
  random_effects <- ranef(mod_fit, condVar = TRUE)[[1]][2]
  slopes <- random_effects[[names(random_effects)[1]]][[hab_var]]
  intercepts <- random_effects[[names(random_effects)[1]]][["(Intercept)"]]        
  
  rnd_eff.df <- data.frame(cbind(intercepts, slopes))
  colnames(rnd_eff.df) <- c("intercept", "slope")
  rnd_eff.df$ID = rownames(random_effects[[1]])
  rownames(rnd_eff.df) <- NULL
  
  return(rnd_eff.df)
}



repeat_func <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  var_across <- vc[vc$grp == "ID", "vcov"]
  var_within <- vc[vc$grp == "ID_season", "vcov"] + var_across
  var_total <- sum(vc$vcov)
  
  rpt_within <- var_within / var_total
  rpt_across <- var_across / var_total
  
  return(c(rpt_within = rpt_within, rpt_across = rpt_across))
}


repeat_func_waBI <- function(fit) {
  
  vc <- as.data.frame(VarCorr(fit))
  var_within <- vc[vc$grp == "ID_season", "vcov"]
  var_total <- sum(vc$vcov)
  
  rpt_within <- var_within / var_total
  rpt_across <- 0
  
  return(c(rpt_within = rpt_within, rpt_across = rpt_across))
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




