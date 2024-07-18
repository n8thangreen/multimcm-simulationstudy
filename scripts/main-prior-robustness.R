# simulation study for manuscript
# cure fraction prior robustness analysis


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

# simulated input data
load(file = "data/input_data.RData")

# global cure fraction hyperparameters
save(params, file = "data/cure_fraction_hyperparameters.RData")


############
# fit model

stan_out <- list()

# i <- 5   # prior
j <- 1   # data

for (i in seq_along(params$mean)) {
  
  data_params <- scenario_data[j, ]
  input_dat <- mutate(input_data[[j]],
                      tx = 1,           # only a single treatment
                      rate = 10^(-10))  # background hazard
  
  ## quick fix by duplicating inputs
  input_dat <- rbind(input_dat,
                     mutate(input_data[[j]],
                            tx = 2,
                            rate = 10^(-100)))
  
  # duplicate for each treatment
  prior_cure <- list(mu_alpha = rep(params$mean[[i]][[1]], 2),
                     sigma_alpha = rep(params$mean[[i]][[2]], 2),
                     mu_sd_cf = rep(params$sd[[i]][[1]], 2),
                     sigma_sd_cf = rep(params$sd[[i]][[2]], 2))
  
  stan_out[[i]] <-
    bmcm_stan(
      input_data = input_dat,
      formula = "Surv(time=times, event=status) ~ 1",
      cureformula = "~ tx + (1 | endpoint)",
      family_latent = data_params$family_latent_model,
      prior_cure = prior_cure,
      centre_coefs = TRUE,
      bg_model = "bg_fixed",
      bg_varname = "rate",
      bg_hr = 1,
      t_max = 5,
      save_stan_code = TRUE)
}

save(stan_out, file = "data/stan_out_cf_priors.RData")

#######
# plot

for (i in seq_along(stan_out)) {
  gg <- plot_S_joint(stan_out[[i]], add_km = TRUE) + xlim(0,5) + facet_wrap(vars(endpoint))
  # ggsave(plot = gg, device = "png", filename = glue::glue("plots/survival_plots_{i}.png"))
}

gg


