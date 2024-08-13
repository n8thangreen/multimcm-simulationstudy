# main_cmdstanr
#
# script to run analysis using cmdstanr

library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)


# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

i <- 1
data <- scenario_data[i, ]
latent_params_true <- eval(parse(text = data$latent_params_true))
n_sim <- 2

sim_params <-
  list(
    n = data$nsample,
    n_endpoints = data$n_endpoints,
    t_cutpoint = data$t_cutpoint,
    mu_cf = data$cf_true,
    sigma_cf = data$sigma_true,
    cf_sample_method = "random",
    cf_indiv = "random",
    distn = data$family_latent_true,
    prop_cens = data$prop_censoring,
    params = latent_params_true)

bmcm_params <- 
  list(
    formula = "Surv(time=times, event=status) ~ 1",
    cureformula = "~ tx + (1 | endpoint)",
    family_latent = data$family_latent_model,
    # duplicate for each treatment
    prior_cure =   
      list(mu_alpha = rep(data$mu_cf_prior, 2),
           sigma_alpha = rep(data$sigma_cf_prior, 2),
           mu_sd_cf = rep(data$mu_sd_cf_prior, 2),
           sigma_sd_cf = rep(data$sigma_sd_cf_prior, 2)),
    centre_coefs = TRUE,
    bg_model = "bg_fixed",
    bg_varname = "rate",
    bg_hr = 1,
    t_max = 5,
    use_cmdstanr = TRUE)

# example simulated data
# for structure
dummy_sample <- do.call(rsurv_cf, sim_params)
dummy_sample$rate <- 10^(-10)
dummy_sample <- 
  mutate(dummy_sample, tx = 1) |> 
  rbind(mutate(dummy_sample, tx = 2))

stan_model <- 
  precompile_bmcm_model(
    input_data = dummy_sample,
    family_latent = bmcm_params$family_latent,
    cureformula = bmcm_params$cureformula,
    use_cmdstanr = TRUE,
    file_path = here::here("stan"))

bmcm_params$precompiled_model_path <- 
  here::here("stan", glue::glue("{stan_model$model_name()}.exe"))

run_scenario(1, sim_params, bmcm_params)

lapply(1:n_sim, \(x) run_scenario(x, sim_params, bmcm_params))




