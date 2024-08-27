# find scenarios where the separate model is different to the hierarchical model
# then use this in the simulation study


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)


# scenario_data
data <- data.frame(
  nsample = 100,
  n_endpoints = 3,
  t_cutpoint = 5,
  mu_cf_prior = -1.39,
  sigma_cf_prior = 0.2,
  mu_sd_cf_prior = 0.4,
  sigma_sd_cf_prior = 2.5,
  cf_true = -1.39,
  sigma_true = 0.4,
  family_latent_true = "weibull",
  family_latent_model = "weibull",
  prop_censoring = 0.5,
  latent_params_true = "list(list(shape = 1, scale = 1), list(shape = 1, scale = 0.1))")

latent_params_true <- eval(parse(text = data$latent_params_true))

prior_cure_hier <- 
  list(mu_alpha = rep(data$mu_cf_prior, 2),
       sigma_alpha = rep(data$sigma_cf_prior, 2),
       mu_sd_cf = rep(data$mu_sd_cf_prior, 2),
       sigma_sd_cf = rep(data$sigma_sd_cf_prior, 2))

# prior parameters for each cure fraction
# duplicate for each treatment
total_sigma <- data$sigma_cf_prior + data$mu_sd_cf_prior
prior_cure_sep <- 
  rep(list(mu_alpha = rep(data$mu_cf_prior, 2),
           sigma_alpha = rep(total_sigma, 2)), data$n_endpoints)

names(prior_cure_sep) <-
  paste0(names(prior_cure_sep), "_", rep(1:data$n_endpoints, each = 2))

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

bmcm_params_hier <- 
  list(
    formula = "Surv(time=times, event=status) ~ 1",
    cureformula = "~ tx + (1 | endpoint)",
    family_latent = data$family_latent_model,
    prior_cure = prior_cure_hier,
    centre_coefs = TRUE,
    bg_model = "bg_fixed",
    bg_varname = "rate",
    bg_hr = 1,
    t_max = 5,
    use_cmdstanr = TRUE)

bmcm_params_sep <- 
  list(
    formula = "Surv(time=times, event=status) ~ 1",
    cureformula = "~ tx + endpoint",
    family_latent = data$family_latent_model,
    prior_cure = prior_cure_sep,
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

# precompiled models

stan_model_hier <- 
  precompile_bmcm_model(
    input_data = dummy_sample,
    family_latent = bmcm_params_hier$family_latent,
    cureformula = bmcm_params_hier$cureformula,
    use_cmdstanr = TRUE,
    file_path = here::here("stan/hierarchical"))

stan_model_sep <- 
  precompile_bmcm_model(
    input_data = dummy_sample,
    family_latent = bmcm_params_sep$family_latent,
    cureformula = bmcm_params_sep$cureformula,
    use_cmdstanr = TRUE,
    file_path = here::here("stan/separate"))

bmcm_params_hier$precompiled_model_path <- stan_model_hier$exe_file()
bmcm_params_sep$precompiled_model_path <- stan_model_sep$exe_file()

# run simulations

run_scenario(1, sim_params, bmcm_params_hier, dir = "output_data/hierarchical/", rstan_format = TRUE)
run_scenario(1, sim_params, bmcm_params_sep, dir = "output_data/separate/", rstan_format = TRUE)


########
# plots

stan_out_hier <- readRDS(here::here("output_data/hierarchical/samples_1.RDS"))
stan_out_sep <- readRDS(here::here("output_data/separate/samples_1.RDS"))

plot_S_joint(stan_out_hier) + xlim(0,5) + facet_wrap(vars(endpoint))
plot_S_joint(stan_out_sep) + xlim(0,5) + facet_wrap(vars(endpoint))

