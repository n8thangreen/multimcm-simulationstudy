# target statistics for each distribution

## median

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

## median cure fraction models

# Weibull distribution
weibull_median_cf <- function(shape, scale, cf) {
  p_star <- (0.5 - cf) / (1 - cf)
  scale * (-log(p_star))^(1 / shape)
}

# Exponential distribution
exp_median_cf <- function(rate, cf) {
  p_star <- (0.5 - cf) / (1 - cf)
  -log(p_star) / rate
}

# Log-normal distribution
lognormal_median_cf <- function(mu, sigma, cf) {
  p_star <- (0.5 - cf) / (1 - cf)
  exp(mu)
}

## rmst

# restricted mean survival time
exp_rmst <- function (rate, tmax) {
  (1 - exp(- rate * tmax))/rate
}

# restricted mean survival time
weibull_rmst <- function (shape, scale, tmax) {
  scale^(-1/shape) * pgamma(scale*tmax^shape, 1/shape + 1) + tmax*exp(-scale*tmax^shape)
}

# restricted mean survival time
gompertz_rmst <- function (shape, scale, tmax) {
  1/scale * (log(1 + shape/scale *(1 - exp(-scale*tmax))) - shape/scale*(1 - exp(-scale*tmax)))
}

# restricted mean survival time
loglogistic_rmst <- function (scale, shape, tmax) {
  exp(-scale/shape) * inc_beta(exp(scale)*tmax^shape/(1 + exp(scale)*tmax^shape), 1 + 1/shape, 1 - 1/shape) + tmax*1/(1 + exp(scale)*tmax^shape)
}

# restricted mean survival time
lognormal_rmst <- function (mu, sigma, tmax) {
  exp(mu + (sigma^2)/2) * Phi((log(tmax) - mu - sigma^2)/sigma) + tmax*(1 - Phi((log(tmax) - mu)/sigma))
}

## rmst of cure fraction models

# restricted mean survival time
exp_rmst_cf <- function (rate, tmax, cf) {
  cf*tmax + (1 - cf)*exp_rmst(rate, tmax)
}

# restricted mean survival time
weibull_rmst_cf <- function (shape, scale, tmax, cf) {
  cf*tmax + (1 - cf)*weibull_rmst(shape, scale, tmax)
}

# restricted mean survival time
gompertz_rmst_cf <- function (shape, scale, tmax, cf) {
  cf*tmax + (1 - cf)*gompertz_rmst(shape, scale, tmax)
}

# restricted mean survival time
loglogistic_rmst_cf <- function (scale, shape, tmax, cf) {
  cf*tmax + (1 - cf)*loglogistic_rmst(scale, shape, tmax)
}

# restricted mean survival time
lognormal_rmst_cf <- function (mu, sigma, tmax, cf) {
  cf*tmax + (1 - cf)*lognormal_rmst(mu, sigma, tmax)
}
