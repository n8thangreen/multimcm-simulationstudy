# performance measures for
# summary estimates

# see
# Using simulation studies to evaluate statistical methods (2017) Tim Morris et al, Stats in Medicine

# targets/estimates:
# * RMST of separate curves
# * cure fraction

# performance measures:
# * bias
# * empirical SE
# * coverage
# * MSE

library(rstan)

load("data/stan_out.RData")


for (i in seq_along(stan_out)) {
  fit <- stan_out[[i]]
  pm[[i]] <- performance_measures(fit, par_nm, true_value)
}


save(pm, file = "data/performance_measures.RData")


########
# plots


