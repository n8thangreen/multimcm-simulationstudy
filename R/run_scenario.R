
# obtain posterior samples for a single simulated data set
#
run_scenario <- function(x, sim_params, bmcm_params, dir = "") {
  set.seed(1234)
  
  # sample data set
  dat <- do.call(rsurv_cf, sim_params)
  
  # background hazard
  dat$rate <- 10^(-10)
  
  # hack because errors with only one treatment
  dat <- 
    mutate(dat, tx = 1) |> 
    rbind(mutate(dat, tx = 2))
  
  bmcm_input <- c(input_data = list(dat), bmcm_params)
  
  # fit model
  fit <- do.call(bmcm_stan, bmcm_input)
  
  samples <- extract_params(fit$output, pattern = "^(median|rmst|cf)") 
  
  # save to file
  path_name <- paste0(dir, "samples_", x, ".csv")
  write.csv(samples, file = path_name)
  
  return()
}

#
extract_params <- function(fit, pattern) {
  if (inherits(fit, "stanfit")) {
    
    samples <- rstan::extract(fit)
    param_names <- grep(pattern, names(samples), value = TRUE)
    extracted_params <- samples[param_names]
    
  } else if (inherits(fit, "CmdStanMCMC")) {
    
    samples <- fit$draws(format = "df")
    param_names <- grep(pattern, names(samples), value = TRUE)
    extracted_params <- samples[, param_names]
  } else {
    stop("Fit object must be of class 'stanfit' (rstan) or 'CmdStanMCMC' (cmdstanr).")
  }
  
  extracted_params
}

