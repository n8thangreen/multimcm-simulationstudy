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


############
# fit model

stan_out <- list()

# i <- 1

for (i in 1:nrow(scenario_data)) {
  
  params <- scenario_data[i, ]
  input_dat <- mutate(input_data[[i]],
                      tx = 1,           # only a single treatment
                      rate = 10^(-10))  # background hazard
  
  ##TODO; errors with single treatment only
  ##      fix in bmcm_stan()
  ## quick fix by duplicating inputs
  input_dat <- rbind(input_dat,
                     mutate(input_data[[i]],
                            tx = 2,
                            rate = 10^(-100)))
  
  # duplicate for each treatment
  prior_cure <- list(mu_alpha = rep(params$mu_cf_prior, 2),
                     sigma_alpha = rep(params$sigma_cf_prior, 2),
                     mu_sd_cf = rep(params$mu_sd_cf_prior, 2),
                     sigma_sd_cf = rep(params$sigma_sd_cf_prior, 2))
  
  stan_out[[i]] <-
    bmcm_stan(
      input_data = input_dat,
      formula = "Surv(time=times, event=status) ~ 1",
      cureformula = "~ tx + (1 | endpoint)",
      family_latent = params$family_latent_model,
      # prior_latent = NA,   ##TODO: how are these used by the Stan code?
      prior_cure = prior_cure,
      centre_coefs = TRUE,
      bg_model = "bg_fixed",
      bg_varname = "rate",
      bg_hr = 1,
      t_max = 5,
      save_stan_code = TRUE)
}

#######
# plot

for (i in seq_along(stan_out)) {
  gg <- plot_S_joint(stan_out[[i]], add_km = TRUE) + xlim(0,5) + facet_wrap(vars(endpoint))
  ggsave(plot = gg, device = "png", filename = glue::glue("plots/survival_plots_{i}.png"))
}

save(stan_out, file = "data/stan_out.RData")


