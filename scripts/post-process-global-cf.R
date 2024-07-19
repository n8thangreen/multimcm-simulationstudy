# performance measures for
# summary estimates
# global cure fraction

library(rstan)
library(purrr)
library(dplyr)
library(ggplot2)


scenario_data <-
  read.csv(here::here("raw-data/scenarios.csv")) |>
  as_tibble() 

# # data scenarios
# file_append <- ""

# prior scenarios
file_append <- "_cf_priors"

load(glue::glue("data/stan_out{file_append}.RData"))

load("data/input_data.RData")

#######################
# performance measures

true_val_names <- c("cf_true", "sigma_true")

# replicate the first scenario values
true_vals <- scenario_data[1, true_val_names] |> 
  rename(lp_cf_global = cf_true, sd_cf = sigma_true) |> 
  # slice(rep(1:n(), each = 5)) |> 
  as.list()

pm_cf <- list()

for (i in seq_along(stan_out)) {
  fit <- stan_out[[i]]
  pm_cf[[i]] <- list()
  
  # estimates
  for (j in c("lp_cf_global", "sd_cf")) {
    
    pm_cf[[i]][[j]] <- bmcm_performance_measures(fit, par_nm = j, true_vals[[j]], append = FALSE)
  }
}

########
# plots

library(tidyr)
library(bayesplot)

# posteriors

# draw grid of plot using ggplot2 with true value marked by vertical line
# for each cure fraction scenario in stan_out

plot_dat <- map(stan_out, "output")

stan_extract <- map(plot_dat, rstan::extract, pars = "sd_cf")
xx <- map(stan_extract, ~as.matrix(.x[[1]])[,1])
xx <- do.call(cbind, xx)
colnames(xx) <- paste0(1:ncol(xx))

p1 <- mcmc_intervals(xx, 
               prob = 0.8,           # 80% credible intervals
               prob_outer = 0.95) +  # 95% credible intervals
geom_vline(xintercept = true_vals$sd_cf, linetype = "dashed", color = "red") +
  xlab("Standard deviation of global cure fraction") +
  theme_minimal()

stan_extract <- map(plot_dat, rstan::extract, pars = "lp_cf_global")
xx <- map(stan_extract, ~as.matrix(.x[[1]])[,1])
xx <- do.call(cbind, xx)
colnames(xx) <- paste0(1:ncol(xx))
p2 <- mcmc_intervals(xx, 
               prob = 0.8,           # 80% credible intervals
               prob_outer = 0.95) +  # 95% credible intervals
geom_vline(xintercept = true_vals$lp_cf_global, linetype = "dashed", color = "red") +
  xlab("Mean of global cure fraction") +
  theme_minimal()

grid_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(grid_plot, filename = glue::glue("plots/posterior_forest_plot_global_cf.png"),
       height = 20, width = 20, dpi = 640, units = "cm")




