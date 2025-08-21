
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



plot_binnedplot <- function(model) {
  
  require(arm)
  
  myplot <- binnedplot(fitted(model), 
                       residuals(model, type = "response"), 
                       nclass = NULL, 
                       xlab = "Expected Values", 
                       ylab = "Average residual", 
                       main = "Binned residual plot", 
                       cex.pts = 0.8, 
                       col.pts = 1, 
                       col.int = "gray")
  
}




# * MODEL OUTPUTS * -------------------------------------------------------

extract_output.lmm <- function(model, ci_level = 0.95) {
  # Get the critical value for the specified CI level
  z <- qnorm(1 - (1 - ci_level) / 2)
  
  # Extract fixed effects summary
  fe_summary <- summary(model)$coefficients
  
  # Convert to data frame
  fixed_effects_df <- as.data.frame(fe_summary)
  colnames(fixed_effects_df) <- c("Estimate", "Std_Error", "df", "t_value", "p_value")
  
  # Calculate confidence intervals
  fixed_effects_df$CI_lower <- fixed_effects_df$Estimate - z * fixed_effects_df$Std_Error
  fixed_effects_df$CI_upper <- fixed_effects_df$Estimate + z * fixed_effects_df$Std_Error
  
  # Add rownames as a column for term names
  fixed_effects_df$Term <- rownames(fixed_effects_df)
  rownames(fixed_effects_df) <- NULL
  
  # Reorder columns
  fixed_effects_df <- fixed_effects_df[, c("Term", "Estimate", "Std_Error", "df", "t_value", "p_value", "CI_lower", "CI_upper")]
  
  return(fixed_effects_df)
}


extract_output.lm <- function(model, ci_level = 0.95) {
  # Get the critical value for the specified CI level
  z <- qnorm(1 - (1 - ci_level) / 2)
  
  # Extract fixed effects summary
  fe_summary <- summary(model)$coefficients
  
  # Convert to data frame
  fixed_effects_df <- as.data.frame(fe_summary)
  colnames(fixed_effects_df) <- c("Estimate", "Std_Error", "t_value", "p_value")
  
  # Calculate confidence intervals
  fixed_effects_df$CI_lower <- fixed_effects_df$Estimate - z * fixed_effects_df$Std_Error
  fixed_effects_df$CI_upper <- fixed_effects_df$Estimate + z * fixed_effects_df$Std_Error
  
  # Add rownames as a column for term names
  fixed_effects_df$Term <- rownames(fixed_effects_df)
  rownames(fixed_effects_df) <- NULL
  
  # Reorder columns
  fixed_effects_df <- fixed_effects_df[, c("Term", "Estimate", "Std_Error", "t_value", "p_value", "CI_lower", "CI_upper")]
  
  return(fixed_effects_df)
}

logit_to_prob <- function(logit_val) {
  1 / (1 + exp(-logit_val)) }

# Helper function to extract emmeans
# get_emmeans <- function(model, by, extra_info) {
#   emmeans(model, reformulate(by)) %>%
#     as.data.frame() %>%
#     mutate(across(everything(), ~.)) %>%
#     bind_cols(extra_info) %>%
#     relocate(names(extra_info), .before = 1)
# }

# Helper function to extract summary output
# get_summary_out <- function(model, model_name, extra_info, varcor = NULL) {
#   s <- summary(model)
#   coefs <- s$coefficients
#   out <- data.frame(
#     intercept     = coefs[1, 1],
#     intercept_sd  = coefs[1, 2],
#     estimate      = coefs[2, 1],
#     se            = coefs[2, 2],
#     df            = ifelse("df" %in% colnames(coefs), coefs[2, "df"], NA),
#     t_value       = coefs[2, "t value"],
#     p_value       = coefs[2, "Pr(>|t|)"],
#     intercept_var = if (!is.null(varcor)) varcor$vcov[1] else NA,
#     resid_var     = if (!is.null(varcor)) varcor$vcov[nrow(varcor)] else NA
#   ) %>% mutate(model = model_name)
#   bind_cols(extra_info, out)
# }

# * PLOTTING FUNCTIONS * --------------------------------------------------

# Plotting individual used estimates
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
    labs(y = "Probability of use")
  
  return(pred_plot)
  
}



save_plot <- function(plot, filename, width = 5, height = 4) {
  png(file = filename, width = width, height = height, units = "in", res = 600)
  print(plot)
  dev.off() }



# * REPEATABILITY CALCULATIONS * ------------------------------------------

 # n_boot = 1000
 # mydata = rep_data.bba.sub
 # myformula = formula

boot_rpt <- function(n_boot, mydata, myformula, nophase = FALSE) {
  
  require(dplyr); require(lme4)
  
  boot_repeatability.annual <- numeric(n_boot)
  boot_repeatability.within <- numeric(n_boot)
  boot_repeatability.between <- numeric(n_boot)
  
  response_var <- all.vars(as.formula(myformula))[1]
  
  # Get unique ring IDs
  unique_rings <- unique(mydata$ring)
  
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  
  for (i in 1:n_boot) {
    
    #print(paste("Bootstrap iteration:", i))
    setTxtProgressBar(pb, i)
    
    boot_rings <- sample(unique_rings, replace = TRUE)
    
    # Resample data with replacement
    boot_data <- do.call(rbind, lapply(boot_rings, function(r) {
      mydata[mydata$ring == r, ]
    }))
    
    boot_data <- boot_data[!is.na(boot_data[[response_var]]), ]
    
    ## Resample if end up with all one colony, or all one phase (unless nophase = TRUE)
    while(n_distinct(boot_data$colony) < 2 | (!nophase && n_distinct(boot_data$phase) < 2)) {
      
      boot_rings <- sample(unique_rings, replace = TRUE)
      
      boot_data <- do.call(rbind, lapply(boot_rings, function(r) {
        mydata[mydata$ring == r, ]
      }))
      
      boot_data <- boot_data[!is.na(boot_data[[response_var]]), ]

    }
    
    # Refit the model
    boot_model <- suppressMessages(lmer(myformula, data = boot_data))
    
    # Extract variance components for bootstrapped model
    boot_var_comp <- as.data.frame(VarCorr(boot_model))
    boot_ring_var <- boot_var_comp$vcov[boot_var_comp$grp == "ring"]
    boot_season_var <- boot_var_comp$vcov[boot_var_comp$grp == "season"]
    boot_ring_season_var <- boot_var_comp$vcov[boot_var_comp$grp == "ring_season"]
    boot_residual_var <- attr(VarCorr(boot_model), "sc")^2
    
    # Calculate within-season repeatability for this bootstrap sample
    boot_repeatability.annual[i] <- boot_season_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
    boot_repeatability.within[i] <- boot_ring_season_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
    boot_repeatability.between[i] <- boot_ring_var / (boot_ring_var + boot_season_var + boot_ring_season_var + boot_residual_var)
  }
  
  close(pb)
  
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




calc_repeatability <- function(focal_var, additional_var, additional_var2, resid_var) {
  return(focal_var / (focal_var + additional_var + additional_var2 + resid_var))
}



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

# extract_rpt <- function(mymodel) {
#   
#   # Extract variance components for bootstrapped model
#   var_comp <- as.data.frame(VarCorr(mymodel))
#   ring_var <- var_comp$vcov[var_comp$grp == "ring"]
#   season_var <- var_comp$vcov[var_comp$grp == "season"]
#   ring_season_var <- var_comp$vcov[var_comp$grp == "ring_season"]
#   residual_var <- attr(VarCorr(mymodel), "sc")^2
#   
#   # Calculate within-season repeatability for this bootstrap sample
#   repeatabilities <- list(
#     rpt_annual = season_var / (ring_var + season_var + ring_season_var + residual_var),
#     rpt_within = ring_season_var / (ring_var + season_var + ring_season_var + residual_var),
#     rpt_between = ring_var / (ring_var + season_var + ring_season_var + residual_var)
#   )
#   
#   rpt_df <- t(data.frame(repeatabilities))
#   colnames(rpt_df) <- "R"
#   return(rpt_df)
# }


extract_rpt <- function(model) {
  # Get variance components
  var_comp <- as.data.frame(VarCorr(model))
  ring_var <- var_comp$vcov[var_comp$grp == "ring"]
  season_var <- var_comp$vcov[var_comp$grp == "season"]
  ring_season_var <- var_comp$vcov[var_comp$grp == "ring_season"]
  residual_var <- attr(VarCorr(model), "sc")^2
  
  # Calculate repeatability estimates
  rpt_annual <- season_var / (ring_var + season_var + ring_season_var + residual_var)
  rpt_within <- ring_season_var / (ring_var + season_var + ring_season_var + residual_var)
  rpt_between <- ring_var / (ring_var + season_var + ring_season_var + residual_var)
  
  # Return the estimates as a vector
  return(c(rpt_annual, rpt_within, rpt_between))
}


extract_var <- function(mymodel) {
  
  my_vars <- as.data.frame(VarCorr(mymodel))
  annual_var = my_vars$vcov[3]
  within_var = my_vars$vcov[1]
  between_var = my_vars$vcov[2]
  
  my_output <- data.frame(rbind(annual_var, within_var, between_var))
  colnames(my_output) <- "var"
  return(my_output)
  
}




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



# Calculate random intercepts and slopes for individual/individual within season

# mod_fit = mod_run.between
# hab_var = "sst_scaled"

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



# n_perm = 1000
# mydata =rep_bba
# myformula =formula
# rptannual =repeatability_output.bba$est[1]
# rptwithin = repeatability_output.bba$est[2]
# rptbetween =repeatability_output.bba$est[3]


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






# * SPATIAL FUNCTIONS * ---------------------------------------------------

## Function to calculate difference in angular bearing
angle_diff <- function(b1, b2) {
  abs_diff <- abs(b1 - b2) # %% 360
  return(ifelse(abs_diff > 180, 360 - abs_diff, abs_diff))
}



angle_diff_rad <- function(b1, b2) {
  abs_diff <- abs(b1 - b2) %% (2 * pi)  # Ensure values wrap correctly
  return(ifelse(abs_diff > pi, (2 * pi) - abs_diff, abs_diff))  # Ensure values are in [0, pi]
}



convert_to_spdf <- function(trip, projection) {
  
  mytrip <- trip %>% ungroup() %>% select(c(longitude, latitude))
  coordinates(mytrip) <- ~longitude+latitude
  proj4string(mytrip) <- proj.dec
  return(mytrip)
}


bearing_from_col <- function(trip, col_lon, col_lat) {
  
  max_location <- trip %>% filter(dist_col == max(dist_col, na.rm = T)) %>%
    ungroup() %>% select(c(longitude, latitude)) %>% filter(row_number() ==1)
  
  return(bearing(c(col_lon, col_lat),
                 c(max_location$longitude, max_location$latitude)))
  
}


# Compute differences and remove any that occur in different breeding phases
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


# df = trip_bearings
# bearing_col_name = "bearing_col"
# comparison_type = "within"
# output_col_name = "bearing_diffs"

calculate_bearing_diffs_sex <- function(df, 
                                    bearing_col_name, 
                                    comparison_type = c("within", "between"), 
                                    output_col_name) {
  
  comparison_type <- match.arg(comparison_type)
  bearing_col_sym <- rlang::sym(bearing_col_name)
  
  # Select relevant columns
  df <- df %>%
    ungroup() %>%
    select(ring, boutID, season, phase, sex, !!bearing_col_sym)
  
  if (comparison_type == "within") {
    
    out <- df %>%
      group_by(ring) %>%
      filter(n() >= 2) %>%
      summarise(
        sex = first(sex),  # Keep sex info
        trip_comparisons = list(combn(boutID, 2, paste, collapse = " vs ")),
        bearing_diff_tmp = list(combn(!!bearing_col_sym, 2, function(x) angle_diff(x[1], x[2]))),
        .groups = "drop" ) %>%
      unnest(cols = c(trip_comparisons, bearing_diff_tmp)) %>%
      rename(!!output_col_name := bearing_diff_tmp) %>%
      separate(trip_comparisons, into = c("tripID_1", "tripID_2"), sep = " vs ") %>%
      left_join(df, by = c("ring", "tripID_1" = "boutID")) %>%
      rename(season_1 = season, phase_1 = phase, sex_1 = sex.x) %>%
      left_join(df, by = c("ring", "tripID_2" = "boutID")) %>%
      rename(season_2 = season, phase_2 = phase, sex_2 = sex) %>%
      mutate(
        comp_ind = "within_ring",
        comp_season = ifelse(season_1 == season_2, "within_season", "between_season"),
        comp_phase = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")
      ) %>%
      filter(comp_phase == "same_phase", sex_1 == sex_2) %>%  # safety check
      mutate(ring_2 = ring) %>%
      rename(ring_1 = ring) %>%
      select(ring_1, tripID_1, season_1, sex_1,
             ring_2, tripID_2, season_2, sex_2,
             !!output_col_name, comp_ind, comp_season)
    
  } else if (comparison_type == "between") {
    
    df1 <- df %>%
      rename_with(~ paste0(.x, "_1"))
    
    df2 <- df %>%
      rename_with(~ paste0(.x, "_2"))
    
    out <- merge(df1, df2, by = NULL) %>%
      filter(ring_1 != ring_2,
             phase_1 == phase_2,
             sex_1 == sex_2) %>%
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
      select(ring_1, boutID_1, season_1, sex_1,
             ring_2, boutID_2, season_2, sex_2,
             !!output_col_name, comp_ind, comp_season) %>%
      rename(tripID_1 = boutID_1, tripID_2 = boutID_2)
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




calculate_within_bearing_diffs.comparison <- function(df, bearing_col_name, output_col_name = "bearing_diff") {
  
  # Convert string column name to symbol
  bearing_col_sym <- rlang::sym(bearing_col_name)
  
  df %>%
    group_by(ring) %>%
    filter(n() >= 2) %>%  # Ensure at least two trips per individual
    summarise(
      trip_comparisons = list(combn(boutID, 2, paste, collapse = " vs ")),
      bearing_diff_tmp = list(combn(!!bearing_col_sym, 2, function(x) angle_diff(x[1], x[2]))),
      .groups = "drop"
    ) %>%
    unnest(cols = c(trip_comparisons, bearing_diff_tmp)) %>%
    rename(!!output_col_name := bearing_diff_tmp) %>%
    separate(trip_comparisons, into = c("tripID_1", "tripID_2"), sep = " vs ") %>%
    left_join(df %>% select(ring, boutID, season, phase), 
              by = c("ring", "tripID_1" = "boutID")) %>%
    rename(season_1 = season, phase_1 = phase) %>%
    left_join(df %>% select(ring, boutID, season, phase), 
              by = c("ring", "tripID_2" = "boutID")) %>%
    rename(season_2 = season, phase_2 = phase) %>%
    mutate(
      comp_ind = "within_ring",
      comp_season = ifelse(season_1 == season_2, "within_season", "between_season"),
      comp_phase = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")
    ) %>%
    arrange(ring, season_1, season_2) %>%
    filter(comp_phase == "same_phase") %>%
    select(-c(phase_1, phase_2, comp_phase)) %>%
    mutate(ring_2 = ring) %>%
    rename(ring_1 = ring) %>%
    relocate(ring_1, tripID_1, season_1, ring_2, tripID_2, season_2,
             !!output_col_name, comp_ind, comp_season) %>%
    mutate(
      tripID_1.tmp = paste0(ring_1, tripID_1),
      tripID_2.tmp = paste0(ring_2, tripID_2),
      comp_ID = paste0(pmin(tripID_1.tmp, tripID_2.tmp), "_", pmax(tripID_1.tmp, tripID_2.tmp))
    ) %>%
    distinct(comp_ID, .keep_all = TRUE) %>%
    select(-c(comp_ID, tripID_1.tmp, tripID_2.tmp))
}




calculate_between_bearing_diffs.comparison <- function(df, bearing_cols = c("bearing_col", "bearing_OG", "bearing_for")) {
  
  # Expand all combinations of IDs and extract components
  pair_df <- expand.grid(ID1 = unique(df$ID),
                         ID2 = unique(df$ID)) %>%
    separate(ID1, into = c("ring_1", "season_1", "dummy1", "tripID_1"), sep = "_", remove = FALSE) %>%
    separate(ID2, into = c("ring_2", "season_2", "dummy2", "tripID_2"), sep = "_", remove = FALSE) %>%
    select(-dummy1, -dummy2) %>%
    mutate(
      tripID_1 = paste0(season_1, "_trip_", tripID_1),
      tripID_2 = paste0(season_2, "_trip_", tripID_2),
      tripID_1.tmp = paste0(ring_1, tripID_1),
      tripID_2.tmp = paste0(ring_2, tripID_2),
      comp_ID = paste0(pmin(tripID_1.tmp, tripID_2.tmp), "_", pmax(tripID_1.tmp, tripID_2.tmp))
    ) %>%
    distinct(comp_ID, .keep_all = TRUE) %>%
    filter(ring_1 != ring_2) %>%
    select(-comp_ID, -tripID_1.tmp, -tripID_2.tmp)
  
  # Start by joining common metadata and phase for filtering
  pair_df <- pair_df %>%
    left_join(df %>% select(ring, boutID, phase), 
              by = c("ring_1" = "ring", "tripID_1" = "boutID")) %>%
    rename(phase_1 = phase) %>%
    left_join(df %>% select(ring, boutID, phase), 
              by = c("ring_2" = "ring", "tripID_2" = "boutID")) %>%
    rename(phase_2 = phase) %>%
    mutate(phase_comp = ifelse(phase_1 == phase_2, "same_phase", "diff_phase")) %>%
    filter(phase_comp == "same_phase") %>%
    select(-phase_1, -phase_2, -phase_comp)
  
  # For each bearing column, compute angle_diff and join
  for (col in bearing_cols) {
    bearing_1 <- paste0(col, "_1")
    bearing_2 <- paste0(col, "_2")
    diff_col  <- paste0(col, ".diff")
    
    tmp <- pair_df %>%
      left_join(df %>% select(ring, boutID, !!sym(col)) %>%
                  rename(ring_1 = ring, tripID_1 = boutID, !!bearing_1 := !!sym(col)),
                by = c("ring_1", "tripID_1")) %>%
      left_join(df %>% select(ring, boutID, !!sym(col)) %>%
                  rename(ring_2 = ring, tripID_2 = boutID, !!bearing_2 := !!sym(col)),
                by = c("ring_2", "tripID_2")) %>%
      mutate(!!diff_col := angle_diff(.data[[bearing_1]], .data[[bearing_2]])) %>%
      select(ring_1, season_1, tripID_1, ring_2, season_2, tripID_2, !!diff_col)
    
    # Join the new difference column to the main df
    pair_df <- left_join(pair_df, tmp, 
                         by = c("ring_1", "season_1", "tripID_1", 
                                "ring_2", "season_2", "tripID_2"))
  }
  
  return(pair_df)
}



deg2rad <- function(deg) {(deg * pi) / (180)}



rad2deg <- function(rad) {(rad * 180) / (pi)}


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



# Vectorized random sampling function
sample_points <- function(geometry, n) {
  st_sample(geometry, size = n, type = "random")
}



# Function to simulate movement from R blog
simulate_movement <- function (sp, env, n, sigma, theta_x, alpha_x, theta_y, alpha_y) {  
  track <- data.frame()  
  track[1,1] <- sp@x  
  track[1,2] <- sp@y  
  for (step in 2:n) {  
    neig <- adjacent(env,   
                     cellFromXY(env, matrix(c(track[step-1,1],  
                                              track[step-1,2]), 1,2)),   
                     directions=8, pairs=FALSE )  
    options <- data.frame()  
    for (i in 1:length(neig)){  
      options[i,1]<-neig[i]  
      options[i,2]<- sp@opt - env[neig[i]]  
    }  
    option <- c(options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ],   
                options[abs(na.omit(options$V2)) == min(abs(na.omit(options$V2))), 1 ])  
    new_cell <- sample(option,1)  
    new_coords <- xyFromCell(env,new_cell)  
    lon_candidate<--9999  
    lat_candidate<--9999  
    
    while ( is.na(extract(env, matrix(c(lon_candidate,lat_candidate),1,2)))) {  
      lon_candidate <- new_coords[1]+ (sigma * rnorm(1)) + (alpha_x * ( theta_x - new_coords[1]))  
      lat_candidate <- new_coords[2]+ (sigma * rnorm(1)) + (alpha_y * ( theta_y - new_coords[2]))  
    }  
    track[step,1] <- lon_candidate  
    track[step,2] <- lat_candidate 
    track$env_value[step] <- extract(env, cellFromXY(env, c(lon_candidate, lat_candidate)))  # Extract environmental value at new position
  }  
  return(track)  
}  
