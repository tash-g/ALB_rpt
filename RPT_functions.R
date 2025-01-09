

# boot_rpt ----------------------------------------------------------------

boot_rpt <- function(n_boot, mydata, myformula) {
  
  boot_repeatability.annual <- numeric(n_boot)
  boot_repeatability.within <- numeric(n_boot)
  boot_repeatability.between <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    
    # Resample data with replacement
    boot_data <- mydata[sample(nrow(mydata), replace = TRUE), ]
    
    # Refit the model
    boot_model <- suppressMessages(lmer(formula, data = boot_data))
    
    # Extract variance components for bootstrapped model
    boot_var_comp <- as.data.frame(VarCorr(boot_model))
    boot_ring_var <- boot_var_comp$vcov[boot_var_comp$grp == "ring"]
    boot_season_var <- boot_var_comp$vcov[boot_var_comp$grp == "season"]
    boot_ring_season_var <- boot_var_comp$vcov[boot_var_comp$grp == "ring:season"]
    boot_residual_var <- attr(VarCorr(boot_model), "sc")^2
    
    # Calculate within-season repeatability for this bootstrap sample
    boot_repeatability.annual[i] <- boot_season_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
    boot_repeatability.within[i] <- boot_ring_season_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
    boot_repeatability.between[i] <- boot_ring_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
  }
  
  # Calculate 95% confidence intervals
  rpt_cis <- list(
    ci_annual = quantile(boot_repeatability.annual, probs = c(0.025, 0.975)),
    ci_within = quantile(boot_repeatability.within, probs = c(0.025, 0.975)),
    ci_between = quantile(boot_repeatability.between, probs = c(0.025, 0.975))  )
  
  ci_annual_str <- paste(signif(rpt_cis$ci_annual[1], 2), signif(rpt_cis$ci_annual[2], 2), sep = ", ")
  ci_within_str <- paste(signif(rpt_cis$ci_within[1], 2), signif(rpt_cis$ci_within[2], 2), sep = ", ")
  ci_between_str <- paste(signif(rpt_cis$ci_between[1], 2), signif(rpt_cis$ci_between[2], 2), sep = ", ")
  
  rpt_cis.df <- data.frame(ci = c(paste0("[", ci_annual_str, "]"),
                                  paste0("[", ci_within_str, "]"),
                                  paste0("[", ci_between_str, "]")))
  
  # Calculate standard errors for within- and between-season repeatability
  rpt_se <- list(
    se_annual = sd(boot_repeatability.annual),
    se_within = sd(boot_repeatability.within),
    se_between = sd(boot_repeatability.between) )
  
  rpt_se.df <- t(data.frame(rpt_se))
  colnames(rpt_se.df) <- "se"
  
  # Calculate mean repeatabilities
  rpt_est <- list(
    est_annual = mean(boot_repeatability.annual),
    est_within = mean(boot_repeatability.within),
    est_between = mean(boot_repeatability.between) )

  rpt_est.df <- t(data.frame(rpt_est))
  colnames(rpt_est.df) <- "est"
  
  # Output
  return(cbind(rpt_est.df, rpt_se.df, rpt_cis.df))

}


# cut_to_breeding ---------------------------------------------------------

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





# calc_repeatability ------------------------------------------------------

calc_repeatability <- function(focal_var, additional_var, additional_var2, resid_var) {
  return(focal_var / (focal_var + additional_var + additional_var2 + resid_var))
}



# calc_rpt_grouped --------------------------------------------------------

calc_rpt_grouped <- function(mydata, myformula, mygroup, group_levels, myspecies) {
  require(dplyr)
  require(lme4)
  
  subset_1 <- subset(mydata, get(mygroup) == group_levels[1])
  subset_2 <- subset(mydata, get(mygroup) == group_levels[2])
  
  mod1 <- lmer(myformula, subset_1)
  mod2 <- lmer(myformula, subset_2)
  
  rpt_1 <- data.frame(extract_rpt(mod1))
  rpt_1$group <- rownames(rpt_1)
  rpt_1$group <- gsub("rpt_", "", rpt_1$group)
  colnames(rpt_1)[1] <- paste0("R_", group_levels[1])
  rownames(rpt_1) <- NULL
  
  rpt_2 <- data.frame(extract_rpt(mod2))
  rownames(rpt_2) <- NULL
  colnames(rpt_2)[1] <- paste0("R_", group_levels[2])
  
  rpt_total <- cbind(rpt_1, rpt_2)
  rpt_total$species <- myspecies
  rpt_total <- rpt_total[,c(4,2,1,3)]
  rpt_total[3:4] <- lapply(rpt_total[3:4], signif, digits = 2)
  
  return(rpt_total)
  
}




# extract_rpt ---------------------------------------------------

extract_rpt <- function(mymodel) {
  
  # Extract variance components for bootstrapped model
  var_comp <- as.data.frame(VarCorr(mymodel))
  ring_var <- var_comp$vcov[var_comp$grp == "ring"]
  season_var <- var_comp$vcov[var_comp$grp == "season"]
  ring_season_var <- var_comp$vcov[var_comp$grp == "ring:season"]
  residual_var <- attr(VarCorr(mymodel), "sc")^2
  
  # Calculate within-season repeatability for this bootstrap sample
  repeatabilities <- list(
    rpt_annual = season_var / (ring_var + season_var + ring_season_var + residual_var),
    rpt_within = ring_season_var / (ring_var + season_var + ring_season_var + residual_var),
    rpt_between = ring_var / (ring_var + season_var + ring_season_var + residual_var)
  )
  
  rpt_df <- t(data.frame(repeatabilities))
  colnames(rpt_df) <- "R"
  return(rpt_df)
}



# extract_var  ------------------------------------------------------------

extract_var <- function(mymodel) {
  
  my_vars <- as.data.frame(VarCorr(mymodel))
  var.annual = my_vars$vcov[3]
  var.within = my_vars$vcov[1]
  var.between = my_vars$vcov[2]
  
  my_output <- data.frame(rbind(annual_var, within_var, between_var))
  colnames(my_output) <- "var"
  return(my_output)
  
}


# extract_variances.brms -------------------------------------------------------

extract_variances.brms <- function(model) {
  
  var_comp <- VarCorr(model)
  
  # Extract variances for Comp.1 and Comp.2
  vars <- list(
    ring_Comp1 = (var_comp$ring$sd[1])^2,
    ringseason_Comp1 = (var_comp$season$sd[1])^2,
    season_Comp1 = (var_comp$`ring:season`$sd[1])^2,
    residual_Comp1 = (var_comp$residual__$sd[1])^2,
    
    ring_Comp2 = (var_comp$ring$sd[2])^2,
    ringseason_Comp2 = (var_comp$season$sd[2])^2,
    season_Comp2 = (var_comp$`ring:season`$sd[2])^2,
    residual_Comp2 = (var_comp$residual__$sd[2])^2
  )
  
  # Calculate repeatability for each component
  repeatabilities <- list(
    repeatability_ring_Comp1 = calc_repeatability(vars$ring_Comp1, vars$season_Comp1, vars$ringseason_Comp1, vars$residual_Comp1),
    repeatability_ringseason_Comp1 = calc_repeatability(vars$ringseason_Comp1, vars$season_Comp1, vars$ring_Comp1, vars$residual_Comp1),
    
    repeatability_ring_Comp2 = calc_repeatability(vars$ring_Comp2, vars$season_Comp2, vars$ringseason_Comp2, vars$residual_Comp2),
    repeatability_ringseason_Comp2 = calc_repeatability(vars$ringseason_Comp2, vars$season_Comp2, vars$ring_Comp2, vars$residual_Comp2)
  )
  
  return(repeatabilities)
}


# get_rnd_effects -------------------------------------------------------------------------

# Calculate random intercepts and slopes for individual/individual within season

get_rnd_effects <- function(mod_fit, hab_var, level) {
  
  ID_var = ifelse(level == "within", "ID_season", "ID")
  
  random_effects <- ranef(mod_fit, condVar = TRUE)
  slopes <- random_effects[[ID_var]][habitat_var]  
  intercepts <- random_effects[[ID_var]]["(Intercept)"]        
  
  rnd_eff.df <- data.frame(cbind(intercepts, slopes))
  colnames(rnd_eff.df) <- c("intercept", "slope")
  rnd_eff.df$ID = rownames(random_effects[[ID_var]])
  rownames(rnd_eff.df) <- NULL
  rnd_eff.df$level = level
  
  return(rnd_eff.df)
}

# loadRdata -----------------------------------------------
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# lrt_test ----------------------------------------------------------------

lrt_test <- function(myformula, rndmeff, mycov, mydata, myspecies) {
  
  # Main model
  mod.main <- lmer(myformula, mydata)
  
  # Null model
  if (rndmeff == "ring") {
    formula_null <- as.formula(paste(mycov, "~ (1|season) + (1|ring:season) + phase + age"))
  } else if (rndmeff == "ring:season") {
    formula_null <- as.formula(paste(mycov, "~ (1|season) + (1|ring) + phase + age"))
  }
  
  ### Fit the null model
  mod.null <- suppressMessages(lmer(formula_null, mydata))
  LRT <- suppressMessages(anova(mod.main, mod.null))
  
  ### Extract results
  LRT_df <- data.frame(LRT)
  
  LRT_df <- LRT_df %>%
    mutate(species = myspecies,
           rndm = ifelse(grepl("null", rownames(LRT_df)), paste0(rndmeff, "_NULL"), rndmeff),
           par = mycov,
           Chisq = round(Chisq, 2),
           Pr..Chisq. = signif(Pr..Chisq., 2)) %>%
    rename(p = Pr..Chisq.) %>%
    select(-c(npar, deviance, Df)) %>%
    relocate(species, par, rndm, .before = AIC)
  rownames(LRT_df) <- NULL
  
  return(LRT_df)
  
}

# perm_repeatability ------------------------------------------------------

perm_repeatability <- function(n_perm, mydata, myformula, rptannual, rptwithin, rptbetween) {
  
  perm_annual <- numeric(n_perm)
  perm_within <- numeric(n_perm)
  perm_between <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    # Shuffle ring labels within each season
    perm_data <- mydata
    perm_data$ring <- ave(perm_data$ring, perm_data$season, FUN = function(x) sample(x))
    
    # Refit the model
    perm_model <- suppressMessages(lmer(formula, data = perm_data))
    
    # Calculate within-season repeatability for this permuted sample
    perm_annual[i] <- extract_rpt(perm_model)[1]
    perm_within[i] <- extract_rpt(perm_model)[2]
    perm_between[i] <- extract_rpt(perm_model)[3]
    
  }
  
  # Calculate the p-value as the proportion of permuted repeatability values
  # greater than or equal to the observed within-season repeatability
  perm_results <- list(
    p_value.annual = mean(perm_annual >= rptannual),
    p_value.within = mean(perm_within >= rptwithin),
    p_value.between = mean(perm_between >= rptbetween)
  )
  
  perm_df <- t(data.frame(perm_results))
  colnames(perm_df) <- "p"
  return(perm_df)
  
}


# plot_ind_predictions ----------------------------------------------------

# Plot individual responses to habitat

plot_ind_predictions <- function(rnd_eff, hab_range, hab_var, mod_fit) {
  
  predictions <- expand.grid(ID = rnd_eff$ID, habitat_val = hab_range)
  
  ## Merge in slope/intercept
  predictions <- merge(predictions, rnd_eff, by = "ID", all.x = T)
  
  predictions$logit <- predictions$intercept + predictions$slope * predictions$habitat_val
  predictions$probability <- 1 / (1 + exp(-predictions$logit))  
  
  ## Get global predictions
  global_intercept <- fixef(mod_fit)["(Intercept)"]
  global_slope <- fixef(mod_fit)[[hab_var]]
  
  predictions.global <- data.frame(habitat_val = hab_range)
  predictions.global$logit <- global_intercept + global_slope * predictions.global$habitat_val
  predictions.global$probability <- 1 / (1 + exp(-predictions.global$logit))  
  
  ## Plot
  pred_plot <- ggplot() +
    geom_line(data = predictions, aes(x = habitat_val, y = probability, colour = ID), alpha = 0.6) +
    geom_line(data = predictions.global, aes(x = habitat_val, y = probability), col = "black",
              linetype = "dashed", linewidth = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Sea surface temperature (°C)", 
         y = "Probability of use")
  
  return(pred_plot)
  
}


# sample_points -----------------------------------------------------------

# Vectorized random sampling function
sample_points <- function(geometry, max_range, n) {
  buffer <- st_buffer(geometry, dist = max_range[1])
  st_sample(buffer, size = n, type = "random")
}