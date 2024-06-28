# simulation study for manuscript


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

# data generation

# # test
# dat <- data.frame(times = rexp(n = 1000,1),
#                   status = 1,
#                   group = 1)
# fit <- survfit(Surv(times, status) ~ group, data = dat)
# plot(fit)

data <- scenario_data[2, ]
latent_params_true <- eval(parse(text = data$latent_params_true))

input_data <-
  rsurv_mix(n = 1000, #data$nsample
            n_endpoints = 100,  #data$n_endpoints,
            t_cutpoint = data$t_cutpoint,
            mu_cf = -1.39,  #data$mu_cf,
            sigma_cf = 0.7,  #data$sigma_true,
            distn = data$family_latent_true,
            prop_cens = data$prop_censoring,
            params = rep(list(latent_params_true[[1]]), 100))

########
# plots
# check Kaplan-Meier

fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

  
  
  
  
  
    
############
# fit model

out <-
  bmcm_stan(
    input_data = input_data,
    formula = "Surv(time=time, event=status) ~ 1",
    cureformula = "~ 1 + (1 | endpoint)",   # single treatment
    family_latent = family_latent,
    centre_coefs = TRUE,
    bg_model = "bg_fixed",
    bg_varname = "rate",
    bg_hr = 1,
    t_max = 5)

