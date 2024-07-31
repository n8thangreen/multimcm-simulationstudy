
#' @import posterior
#' 
performance_measures <- function(theta_hat, true_value,
                                 theta_hat_low = NA, theta_hat_upp = NA) {
  # estimates
  bias <- mean(theta_hat, na.rm = TRUE) - true_value
  empirical_se <- stats::sd(theta_hat, na.rm = TRUE)
  mse <- mean((theta_hat - true_value)^2, na.rm = TRUE)
  coverage <- calc_coverage(true_value, theta_hat_low, theta_hat_upp)
  
  # Monte Carlo SE of estimate
  mc_se_bias <- stats::sd(theta_hat, na.rm = TRUE) / sqrt(length(theta_hat))
  mc_se_empirical_se <- stats::sd(theta_hat, na.rm = TRUE) / sqrt(2 * (length(theta_hat) - 1))
  mc_se_coverage <- sqrt(coverage * (1 - coverage) / length(theta_hat))
  # mc_se_mse ##TODO 
  
  c(bias = bias,
    empirical_se = empirical_se,
    mse = mse,
    coverage = coverage)
}

#' @param theta.hat.low lower bound of interval estimate
#' @param theta.hat.upp upper bound of interval estimate
calc_coverage <- function(true_value,
                          theta_hat_low = NA, theta_hat_upp = NA) {
  if (any(is.na(theta_hat_low))) return()
  
  nsim <- length(theta_hat_low)
  bin_vals <- ifelse(true_value >= theta_hat_low & true_value <= theta_hat_upp, 1, 0)
  
  sum(bin_vals)/nsim
}

#' statistics for all endpoints
#' for deterministic cure fraction data
#' 
#' @importFrom posterior merge_chains as_draws
#' 
bmcm_performance_measures <- function(fit, par_nm, true_vals, append = TRUE) {
  res <- NULL
  n_endpoints <- length(true_vals)
  
  for (i in seq_len(n_endpoints)) {
    if (append) {
      par_nm_ <- paste0(par_nm, "_", i)
    } else {
      par_nm_ <- par_nm
    }
    true_val <- true_vals[i]
    
    # extract posterior samples
    stan_extract <- rstan::extract(fit$output, pars = par_nm_)
    all_samples <- merge_chains(as_draws(stan_extract))[[1]]
    samples <- all_samples[[par_nm_]]
    
    res <- rbind(res, performance_measures(samples, true_val))
  }
  
  res
}

#' statistics for all endpoints
#' for full probabilistic analysis
#' 
#' @importFrom posterior merge_chains as_draws
#' 
bmcm_performance_measures_N <- function(stan_out_list, par_nm, true_vals, append = TRUE) {

  res <- NULL
  n_endpoints <- length(true_vals)
  
  for (i in seq_len(n_endpoints)) {
    if (append) {
      par_nm_ <- paste0(par_nm, "_", i)
    } else {
      par_nm_ <- par_nm
    }
    true_val <- true_vals[i]
    
    # extract posterior samples
    samples <- lapply(stan_out_list,
                      function(x) {
                        stan_extract <- rstan::extract(x$output, pars = par_nm_)
                        all_samples <- merge_chains(as_draws(stan_extract))[[1]]
                        all_samples[[par_nm_]]
                      })
    theta_hat <- sapply(samples, mean)
    theta_hat_low <- sapply(samples, quantile, probs = 0.025, na.rm = TRUE)
    theta_hat_upp <- sapply(samples, quantile, probs = 0.975, na.rm = TRUE)

    ##TODO: deal with NaN?    
    # theta_hat_low <- sapply(samples, function(x) {
    #   if (all(is.nan(x))) {
    #     NaN
    #   } else {
    #     quantile(x, probs = 0.025, na.rm = TRUE)
    #   }
    # })
    
    res <- rbind(res,
                 performance_measures(
                   theta_hat, true_val, theta_hat_low, theta_hat_upp))
  }
  
  res
}
