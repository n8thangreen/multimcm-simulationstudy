# performance measures for
# summary estimates

# see
# Using simulation studies to evaluate statistical methods (2017) Tim Morris et al, Stats in Medicine

##TODO: 
## * take average across curves?

# targets/estimates:
# * RMST of separate curves
# * cure fraction
# * median survival time 

# performance measures:
# * bias
# * empirical SE
# * coverage
# * MSE


##TODO: how to use posterior package?

library(rstan)

load("data/stan_out.RData")
load("data/input_data.RData")

target_names <- c("median")  #, "cf", "rmst")
pm <- list()

# scenarios
for (i in seq_along(stan_out)) {
  fit <- stan_out[[i]]
  
  # estimates
  for (j in target_names) {
    true_vals <- attr(input_data[[i]], which = j)
    pm[[i]] <- bmcm_performance_measures(fit, par_nm = j, true_vals)
  }
}


save(pm, file = "data/performance_measures.RData")


########
# plots


#########
# tables
