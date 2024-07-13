
#' @import posterior
#' 
performance_measures <- function(samples, true_value, ci_level = 0.95) {
  
  bias <- mean(samples) - true_value
  empirical_se <- sd(samples)
  mse <- mean((samples - true_value)^2)
  
  lower_quantile <- (1 - ci_level) / 2
  upper_quantile <- 1 - lower_quantile
  credible_interval <- quantile(samples, probs = c(lower_quantile, upper_quantile))
  coverage <- ifelse(true_value >= credible_interval[1] &
                       true_value <= credible_interval[2],
                     1, 0)
  
  list(
    bias = bias,
    empirical_se = empirical_se,
    mse = mse,
    coverage = unname(coverage))
}

# median survival time for a Weibull distribution
weibull_median <- function(shape, scale) {
  scale * (log(2))^(1 / shape)
}

# median survival time for an Exponential distribution
exp_median <- function(rate) {
  log(2) / rate
}

# median survival time for a Log-normal distribution
lognormal_median <- function(mu, sigma) {
  exp(mu)
}


# statistics for all endpoints
bmcm_performance_measures <- function(fit, par_nm, true_vals) {
  
  res <- list()
  n_endpoints <- fit$formula$cure$cf_idx
  
  for (i in seq_len(n_endpoints)) {
    par_nm_ <- paste0(par_nm, "_", i)
    true_val <- true_vals[i]
    
    # extract posterior samples
    stan_extract <- rstan::extract(fit$output, pars = par_nm)
    all_samples <- merge_chains(as_draws(stan_extract))[[1]]
    samples <- all_samples[[par_nm]]
    
    res[[i]] <- performance_measures(samples, true_val)
  }
  
  res
}
