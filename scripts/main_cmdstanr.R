# main_cmdstanr
#
# script to run analysis using cmdstanr

library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)


#
run_scenario <- function(x, sim_params, bmcm_params) {
  
  # sample data set
  dat <- do.call(rsurv_cf, sim_params)
  
  # fit model
  fit <- do.call(bmcm_stan, list(dat, bmcm_params))
  
  # extract posterior samples
  samples <- fit$df_draws
  
  # save to file
  path_name <- paste0("samples_", x, ".csv")
  write.csv(samples, file = path_name)
  
  invisible(samples)
}


# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

i <- 1
data <- scenario_data[i, ]
latent_params_true <- eval(parse(text = data$latent_params_true))
n_sim <- 10

params <- list(n = data$nsample,
               n_endpoints = data$n_endpoints,
               t_cutpoint = data$t_cutpoint,
               mu_cf = data$cf_true,
               sigma_cf = data$sigma_true,
               cf_sample_method = "random",
               cf_indiv = "random",
               distn = data$family_latent_true,
               prop_cens = data$prop_censoring,
               params = latent_params_true)

lapply(1:n_sim, \(x) run_scenario(x, params))




