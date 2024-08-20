
#' @import posterior
#' 
performance_measures <- function(theta_hat, true_value,
                                 theta_hat_low = NA, theta_hat_upp = NA) {
  # estimates
  bias <- mean(theta_hat - true_value, na.rm = TRUE)
  relative_bias <- mean((theta_hat - true_value)/true_value, na.rm = TRUE)
  empirical_se <- stats::sd(theta_hat, na.rm = TRUE)
  mse <- mean((theta_hat - true_value)^2, na.rm = TRUE)
  coverage <- calc_coverage(true_value, theta_hat_low, theta_hat_upp)
  
  # Monte Carlo SE of estimate
  mc_se_bias <- stats::sd(theta_hat, na.rm = TRUE) / sqrt(length(theta_hat))
  mc_se_empirical_se <- stats::sd(theta_hat, na.rm = TRUE) / sqrt(2 * (length(theta_hat) - 1))
  mc_se_coverage <- sqrt(coverage * (1 - coverage) / length(theta_hat))
  # mc_se_mse ##TODO 
  
  c(bias = bias,
    relative_bias = relative_bias,
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
#' computed locally
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
    
    # clean means
    which_keep <- !(theta_hat > 10 | is.infinite(theta_hat) | is.na(theta_hat))
    
    theta_hat <- theta_hat[which_keep]
    theta_hat_low <- theta_hat_low[which_keep]
    theta_hat_upp <- theta_hat_upp[which_keep]
    true_val <- true_val[which_keep]
    
    res <- rbind(res,
                 performance_measures(
                   theta_hat, true_val, theta_hat_low, theta_hat_upp))
  }
  
  res
}


#' Performance measures from cluster
#' 
#' statistics for all endpoints for full probabilistic analysis
#' using posterior samples obtain using cluster
#' 
#' @param stan_out_list 
#' @param par_nm 
#' @param true_vals 
#' @param append when there are multiple curves near to add a number after par_nm
#'
#' @importFrom posterior merge_chains as_draws
#' 
performance_measures_cluster <- function(stan_out_list,
                                         par_nm,
                                         true_vals,
                                         append = TRUE) {
  res <- NULL
  n_endpoints <- nrow(true_vals[[1]])
  
  for (i in seq_len(n_endpoints)) {
    
    theta_stats <- samples_summary_stats(
      stan_out_list, par_nm, true_vals, i, append)
    
    res <- rbind(res,
                 performance_measures(
                   theta_stats$theta_hat,
                   theta_stats$theta_true,
                   theta_stats$theta_hat_low,
                   theta_stats$theta_hat_upp))
  }
  
  res
}

#
samples_summary_stats <- function(stan_out_list,
                                  par_nm,
                                  true_vals = NA,
                                  endpoint_id,
                                  append = TRUE) {
  if (append) {
    par_nm_ <- paste0(par_nm, "_", endpoint_id)
  } else {
    par_nm_ <- par_nm
  }
  
  if (any(!is.na(true_vals))) {
    theta_true <- sapply(true_vals, \(x) x[endpoint_id, par_nm])
  } else {
    theta_true <- NA
  }
  
  # extract posterior samples
  samples <- purrr::map(stan_out_list, par_nm_)
  
  theta_hat <- sapply(samples, mean)
  theta_hat_low <- sapply(samples, quantile, probs = 0.025, na.rm = TRUE)
  theta_hat_upp <- sapply(samples, quantile, probs = 0.975, na.rm = TRUE)
  
  # clean means
  which_keep <- !(theta_hat > 10 | is.infinite(theta_hat) | is.na(theta_hat))
  
  theta_hat <- theta_hat[which_keep]
  theta_hat_low <- theta_hat_low[which_keep]
  theta_hat_upp <- theta_hat_upp[which_keep]
  theta_true <- theta_true[which_keep]
  
  data.frame(
    theta_hat,
    theta_hat_low,
    theta_hat_upp,
    theta_true)
}

