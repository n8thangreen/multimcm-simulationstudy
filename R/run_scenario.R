
#' Run scenario
#' 
#' obtain posterior samples for a single simulated data set
#' @param x sample number
#' @param rstan_format list (rstan) or flat dataframe (cmdstanr) format; default FALSE
#' @param ... additional arguments passed to \code{\link{bmcm_stan}}
#'
run_scenario <- function(x, sim_params, bmcm_params,
                         dir = "",
                         seed = NULL,
                         rstan_format = FALSE,
                         ...) {
  set.seed(seed)
  dots <- list(...)
  
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
  
  bmcm_input <- c(input_data = list(dat), bmcm_params, dots)
  
  # fit model
  fit <- do.call(bmcm_stan, bmcm_input)
  
  # samples <- extract_params(fit$output, pattern = "^(median|rmst|cf)", rstan_format)
  ## for testing
  samples <- extract_params(fit$output, pattern = "^(median|rmst|cf|sd_cf)", rstan_format)
  
  # save to file
  if (rstan_format) {
    path_name <- paste0(dir, "samples_", x, ".RDS")
    saveRDS(samples, file = path_name)
    
    input_data_path <- paste0(dir, "true_values_", x, ".RDS")
    saveRDS(true_values, file = input_data_path)
  } else {
    path_name <- paste0(dir, "samples_", x, ".csv")
    write.csv(samples, file = path_name)
    
    input_data_path <- paste0(dir, "true_values_", x, ".csv")
    write.csv(true_values, file = input_data_path)
  }
  
  return()
}

#' @param fit Stan model object
#' 
extract_params <- function(fit, pattern = "", rstan_format = FALSE) {
  if (inherits(fit, "stanfit")) {
    
    samples <- rstan::extract(fit)
    param_names <- grep(pattern, names(samples), value = TRUE)
    extracted_params <- samples[param_names]
    
  } else if (inherits(fit, "CmdStanMCMC")) {
    
    if (rstan_format) {
      # list format
      samples <- as_draws_rvars(fit$draws())
      samples_df <- lapply(samples, \(x) as_draws_df(x))
      
      # remove extra columns
      samples_df <-
        lapply(samples_df, function(x) {
          drop_cols <- colnames(x) %in% c(".chain", ".iteration", ".draw")
          x[, !drop_cols]
        })
      
      param_names <- grep(pattern, names(samples_df), value = TRUE)
      extracted_params <- samples_df[param_names]
    } else {
      # flat rectangular data frame
      samples <- fit$draws(format = "df")
      param_names <- grep(pattern, names(samples), value = TRUE)
      extracted_params <- samples[, param_names]
    }
  } else {
    stop("Fit object must be of class 'stanfit' (rstan) or 'CmdStanMCMC' (cmdstanr).")
  }
  
  extracted_params
}

