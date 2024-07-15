
#' @import posterior
#' 
performance_measures <- function(samples, true_value, ci_level = 0.95) {
  
  bias <- mean(samples, na.rm = TRUE) - true_value
  empirical_se <- stats::sd(samples, na.rm = TRUE)
  mse <- mean((samples - true_value)^2, na.rm = TRUE)
  coverage <- calc_coverage(samples, true_value, ci_level)
  
  c(bias = bias,
    empirical_se = empirical_se,
    mse = mse,
    coverage = coverage)
}

#
calc_coverage <- function(samples, true_value, ci_level = 0.95) {
  
  lower_quantile <- (1 - ci_level) / 2
  upper_quantile <- 1 - lower_quantile
  credible_interval <- quantile(samples, probs = c(lower_quantile, upper_quantile), na.rm = TRUE)
  coverage <- ifelse(true_value >= credible_interval[1] &
                       true_value <= credible_interval[2],
                     1, 0)
  unname(coverage)
}

#' statistics for all endpoints
#' 
#' @importFrom posterior merge_chains as_draws
#' 
bmcm_performance_measures <- function(fit, par_nm, true_vals) {
  res <- NULL
  n_endpoints <- length(true_vals)
  
  for (i in seq_len(n_endpoints)) {
    par_nm_ <- paste0(par_nm, "_", i)
    true_val <- true_vals[i]
    
    # extract posterior samples
    stan_extract <- rstan::extract(fit$output, pars = par_nm_)
    all_samples <- merge_chains(as_draws(stan_extract))[[1]]
    samples <- all_samples[[par_nm_]]
    
    res <- rbind(res, performance_measures(samples, true_val))
  }
  
  res
}
