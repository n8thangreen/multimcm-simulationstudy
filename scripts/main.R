# simulation study for manuscript


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

##################
# data generation

# # elicit latent survival curves
# dat <- data.frame(times = rexp(n = 1000,1),
#                   status = 1,
#                   group = 1)
# fit <- survfit(Surv(times, status) ~ group, data = dat)
# plot(fit)

# # elicit sd of cure fraction prior
# # on the logistic scale
#
# sigma_cf <- 1     # uninformative
# sigma_cf <- 0.2   # informative
# alpha <- rnorm(n = 10000, mean = -1.39, sigma_cf)
# 
# hist(exp(alpha)/(1 + exp(alpha)), breaks = 50, xlim = c(0,1))


input_data <- list()

for (i in 1:16) {
  
  data <- scenario_data[i, ]
  latent_params_true <- eval(parse(text = data$latent_params_true))
  
  input_data[[i]] <-
    rsurv_mix(n = data$nsample,
              n_endpoints = data$n_endpoints,
              t_cutpoint = data$t_cutpoint,
              mu_cf = data$cf_true,
              sigma_cf = data$sigma_true,
              cf_sample_method = "quantiles",
              distn = data$family_latent_true,
              prop_cens = data$prop_censoring,
              params = latent_params_true)
}

save(input_data, file = "data/input_data.RData")
# load(file = "data/input_data.RData")


########
# plots

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data[[1]])

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

png(filename = "plots/simulated_survival_plots.png", width = 20, height = 20, units = "cm", res = 300)

# grid of plots
par(mfrow = c(4,4))
# par(mfrow = c(1,1))
fit_dat <- list()

for (i in seq_along(input_data)) {
  fit_dat[[i]] <- survfit(Surv(times, status) ~ endpoint, data = input_data[[i]])
  
  data <- scenario_data[i, ]
  plot_title <- glue("{data$data_id}: e={data$n_endpoints} n={data$nsample} p={data$prop_censoring} s={data$sigma_true}")

  cure_fractions <-
    input_data[[i]] |>
    group_by(endpoint) |>
    summarise(x = sum(curestatus == 2) / n())
  
  plot(fit_dat[[i]], col = rgb(0,0,0, alpha = 0.3),
       main = plot_title,
       xlab = "Time", ylab = "Survival probability")
  abline(h = cure_fractions$x, col = "pink")
  abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")
}
dev.off()

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


