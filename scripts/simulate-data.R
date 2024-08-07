# simulation study for manuscript
#
# generate data for deterministic and probabilistic
# cure fraction scenarios


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

n_scenarios <- 16

## deterministic cure fraction sampling

for (i in 1:n_scenarios) {
  data <- scenario_data[i, ]
  latent_params_true <- eval(parse(text = data$latent_params_true))
  
  input_data[[i]] <-
    rsurv_cf(n = data$nsample,
             n_endpoints = data$n_endpoints,
             t_cutpoint = data$t_cutpoint,
             mu_cf = data$cf_true,
             sigma_cf = data$sigma_true,
             cf_sample_method = "quantiles",
             cf_indiv = "fixed",
             distn = data$family_latent_true,
             prop_cens = data$prop_censoring,
             params = latent_params_true)
}

save(input_data, file = "data/determ_input_data.RData")
# load(file = "data/determ_input_data.RData")

## full probabilistic simulation

input_data <- list()
n_sim <- 10

for (i in 1:n_scenarios) {
  print(i)
  
  data <- scenario_data[i, ]
  latent_params_true <- eval(parse(text = data$latent_params_true))
  
  rsurv_args <- list(n = data$nsample,
                     n_endpoints = data$n_endpoints,
                     t_cutpoint = data$t_cutpoint,
                     mu_cf = data$cf_true,
                     sigma_cf = data$sigma_true,
                     cf_sample_method = "random",
                     cf_indiv = "random",
                     distn = data$family_latent_true,
                     prop_cens = data$prop_censoring,
                     params = latent_params_true)
  
  for (j in 1:n_sim) {
    suppressMessages(
      input_data[[j]] <- do.call(rsurv_cf, rsurv_args)
    )
  }
  
  # save each scenario to separate file
  # good idea for large n_sim
  file_name <- glue::glue("N{data$nsample}_ne{data$n_endpoints}_pcens{data$prop_censoring}_sigma{data$sigma_true}")
  save(input_data, file = paste0("data/", file_name, ".RData"))
}


########
# plots

## deterministic cure fraction

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data[[1]])

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

png(filename = "plots/simulated_survival_plots.png", width = 20, height = 20, units = "cm", res = 300)

# grid of plots
x11()
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

