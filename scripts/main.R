# simulation study for manuscript


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

# data generation

# test
# dat <- data.frame(times = rexp(n = 1000,1),
#                   status = 1,
#                   group = 1)
# fit <- survfit(Surv(times, status) ~ group, data = dat)
# plot(fit)

data <- scenario_data[1,]
latent_params_true <- eval(parse(text = data$latent_params_true))

input_data <-
  rsurv_mix(n = 1000,
            n_endpoints = data$n_endpoints,
            t_cutpoint = data$t_cutpoint,
            mu_cf = data$mu_cf,
            sigma_cf = data$sigma_true,
            distn = data$family_latent_true,
            prop_cens = data$prop_censoring,
            params = latent_params_true)

########
# plots
# check Kaplan-Meier

fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat)

  
  
  
  
  
    
############
# fit model

out <-
  bmcm_stan(
    input_data = input_data,
    formula = "Surv(time=time, event=status) ~ 1",
    cureformula = "~ treat + (1 | endpoint)",
    family_latent = family_latent,
    centre_coefs = TRUE,
    bg_model = "bg_fixed",
    bg_varname = "rate",
    bg_hr = 1,
    t_max = 400)

