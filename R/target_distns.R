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

## rmst

# restricted mean survival time
exp_rmst <- function (rate, tmax) {
  (1 - exp(- rate * tmax))/rate
}

# restricted mean survival time
weibull_rmst <- function (shape, scale, tmax) {
  scale^(-1/shape) * gamma_q(scale*tmax^shape, 1/shape + 1) + tmax*exp(-scale*tmax^shape)
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

