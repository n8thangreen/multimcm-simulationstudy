
#
performance_measures <- function(fit, par_nm, true_value) {
  
  # Extract posterior samples
  samples <- extract(fit, pars = par_nm)$param
  
  bias <- mean(samples) - true_value
  empirical_se <- sd(samples)
  mse <- mean((samples - true_value)^2)
  
  ci_level <- 0.95
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
    coverage = coverage)
}
