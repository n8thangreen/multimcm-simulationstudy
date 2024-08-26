
#' run_scenario
#' 
#' obtain posterior samples for a single simulated data set
#' @param x sample number
#'
run_scenario <- function(x, sim_params, bmcm_params, dir = "", seed = NULL) {
  set.seed(seed)
  
  # sample data set
  dat <- do.call(rsurv_cf, sim_params)
  
  true_values <- data.frame(
    rmst = attr(dat, which = "rmst"),
    median = attr(dat, which = "median"),
    cf = attr(dat, which = "cf"))
  
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

  input_data_path <- paste0(dir, "true_values_", x, ".csv")
  write.csv(true_values, file = input_data_path)
  
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

