# find scenarios where the separate model is different to the hierarchical model
# then use this in the simulation study


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)


############
# functions 

invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

#
truncated_cauchy <- function(n, location, scale, a) {
  samples <- numeric(n)
  i <- 1
  while (i <= n) {
    x <- rcauchy(1, location, scale)
    if (x > a) {
      samples[i] <- x
      i <- i + 1
    }
  }
  samples
}

################
# cure fraction

## true
invlogit(rnorm(1000, -1, 0.1)) |> 
  density() |> plot()

## prior
n <- 10000

# weakly informative
# mu <- rnorm(n, mean = -1, sd = 0.7)
mu <- rnorm(n, mean = 0, sd = 0.7)
sd <- truncated_cauchy(n, location = 0.05, scale = 0.1, 0)  # between-group

# # very informative
# mu <- rnorm(n, mean = -1, sd = 0.01)
# sd <- truncated_cauchy(n, location = 0.05, scale = 0.05, 0)

mu_cf <- invlogit(mu)     # global cure fraction
sd <- sd[sd < 10]

density(mu) |> plot()
density(mu_cf) |> plot()
density(sd) |> plot(xlim = c(0,1))
abline(v = median(sd), col = "red")

# group-level cure fraction
rnorm(n = n, mean = mu, sd = sd) |>
  invlogit() |> 
  density() |>
  plot(lty = 2, col = "red") #, xlim = c(0, 0.6))
density(mu_cf) |> lines()

#########################
# latent survival curves

n <- 10000

# less informative
shape <- rgamma(n, shape = 10, scale = 0.1)
scale1 <- rlnorm(n, meanlog = 0, sdlog = 0.1)
scale4 <- rlnorm(n, meanlog = 1.4, sdlog = 0.02)

# very informative
# shape <- rgamma(n, shape = 1000, scale = 0.001)
# scale1 <- rlnorm(n, meanlog = 0, sdlog = 0.01)         # scale = 1
# scale4 <- rlnorm(n, meanlog = 1.387, sdlog = 0.002)    # scale = 4
# log_scale4 <- rnorm(n, mean = 1.387, sd = 0.002)
# exp_scale4 <- exp(log_scale4)

density(shape) |> plot()
density(scale1) |> plot()
density(scale4) |> plot()
# density(exp_scale4) |> plot()

plot(NULL, ylim = c(0,1), xlim = c(0,5))
for (i in 1:100) {
  pweibull(q = seq(0, 5, 0.1), shape = shape[i], scale = scale1[i], lower.tail = FALSE) |> 
    lines(x = seq(0, 5, 0.1), col = "grey")
}
# true survival curve
pweibull(q = seq(0, 5, 0.1), shape = 1, scale = 1, lower.tail = FALSE) |> 
  lines(x = seq(0, 5, 0.1), col = "red")

plot(NULL, ylim = c(0,1), xlim = c(0,5))
for (i in 1:100) {
  pweibull(q = seq(0,5,0.1), shape = shape[i], scale = scale4[i], lower.tail = FALSE) |> 
    lines(x = seq(0,5,0.1), col = "grey")
}
# true survival curve
pweibull(q = seq(0, 5, 0.1), shape = 1, scale = 4, lower.tail = FALSE) |> 
  lines(x = seq(0, 5, 0.1), col = "red")


###########
# analysis
###########

# scenario_data
data <- data.frame(
  nsample = 100,
  n_endpoints = 10,
  t_cutpoint = 5,
  mu_cf_prior = 0,
  # mu_cf_prior = -1,        # true value
  sigma_cf_prior = 0.01,   # very informative
  # sigma_cf_prior = 0.7,
  mu_sd_cf_prior = 0.05,     # very informative
  # sigma_sd_cf_prior = 0.05,
  sigma_sd_cf_prior = 0.1,   # weaker
  cf_true = -1,
  sigma_true = 0.1,
  family_latent_true = "weibull",
  family_latent_model = "weibull",
  prop_censoring = 0,
  latent_params_true =
    "list(
     list(shape = 1, scale = 4),
     list(shape = 1, scale = 1))",  # scale is 1/rate
  latent_shape_prior = 
    "list(
     list(a = 1000, b = 1/0.01),
     list(a = 1000, b = 1/0.01))",
     # list(a = 1000, b = 1/0.001),   # very informative
     # list(a = 1000, b = 1/0.001))",
  latent_scale_prior =
    "list(
     list(mu = 1.4, sigma = 0.02),
     list(mu = 0, sigma = 0.1))"
     # list(mu = 1.387, sigma = 0.002),  # very informative
     # list(mu = 0, sigma = 0.01))"
)

# convert to list
latent_params_true <- eval(parse(text = data$latent_params_true))
latent_shape_prior <- eval(parse(text = data$latent_shape_prior))
latent_scale_prior <- eval(parse(text = data$latent_scale_prior))

prior_cure_hier <- 
  list(mu_alpha = rep(data$mu_cf_prior, 2),
       sigma_alpha = rep(data$sigma_cf_prior, 2),
       mu_sd_cf = rep(data$mu_sd_cf_prior, 2),
       sigma_sd_cf = rep(data$sigma_sd_cf_prior, 2))

# prior parameters for each cure fraction
# duplicate for each treatment
total_sigma <- sqrt(data$sigma_cf_prior^2 + data$mu_sd_cf_prior^2 + data$sigma_sd_cf_prior^2)

# # check approximately equivalent to hierarchical
# rnorm(n = n, mean = mu, sd = sd) |>
#   density() |> plot()
# rnorm(n = n, mean = -1, sd = total_sigma) |> 
#   density() |> lines(col = "red")

prior_cure_sep <-
  rep(list(mu_alpha = rep(data$mu_cf_prior, 2),
           sigma_alpha = rep(total_sigma, 2)), data$n_endpoints)

# ## for testing
# prior_cure_sep <- 
#   rep(list(mu_alpha = rep(0, 2),
#            sigma_alpha = rep(2, 2)), data$n_endpoints)

names(prior_cure_sep) <-
  paste0(names(prior_cure_sep), "_", rep(1:data$n_endpoints, each = 2))

# need to convert to array for stan
latent_scale_prior <- map(latent_scale_prior, ~ map(.x, as.array))

##TODO: move this to inside bmcm_stan()
# prior parameters for each latent curve
# for each end point
# recycle is length of parameters < n_endpoints

prior_latent_list <- 
  c(list_flatten(rep(latent_shape_prior, length.out = data$n_endpoints)),
    list_flatten(rep(latent_scale_prior, length.out = data$n_endpoints)))


# parameter name
names(prior_latent_list) <-
  paste0(names(prior_latent_list), rep(c("_shape", "_S"), each = 2*data$n_endpoints))

# endpoint index
names(prior_latent_list) <-
  paste0(names(prior_latent_list), "_", rep(1:data$n_endpoints, each = 2))

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
    prior_latent = prior_latent_list,
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
    prior_latent = prior_latent_list,
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

run_scenario(1, sim_params, bmcm_params_hier, dir = "output_data/hierarchical/", rstan_format = FALSE,
             iter_warmup = 200,
             iter_sampling = 1000,
             save_warmup = FALSE,
             thin = 1) #, seed = 1234)

run_scenario(1, sim_params, bmcm_params_sep, dir = "output_data/separate/", rstan_format = FALSE,
             iter_warmup = 200,
             iter_sampling = 1000,
             save_warmup = FALSE,
             thin = 1) #, seed = 1234)

########
# plots
########

# stan_out_hier <- readRDS(here::here("output_data/hierarchical/samples_1.RDS"))
# stan_out_sep <- readRDS(here::here("output_data/separate/samples_1.RDS"))
stan_out_hier <- read.csv(here::here("output_data/hierarchical/samples_1.csv"))
stan_out_sep <- read.csv(here::here("output_data/separate/samples_1.csv"))

stan_out_hier_true <- read.csv(here::here("output_data/hierarchical/true_values_1.csv"))
stan_out_sep_true <- read.csv(here::here("output_data/separate/true_values_1.csv"))

param_names <- names(stan_out_sep)
cf_names <- param_names[-1]  # remove row number
cf_names <- sort(cf_names[!grepl(pattern = "2\\.$", cf_names)])

# cure fraction separate prior
cf_sep_lin <- rnorm(n = 10000, mean = data$mu_cf_prior, sd = total_sigma)
cf_data <- data.frame(cf_sep = invlogit(cf_sep_lin))

plot_list <- list()

for (i in cf_names) {
  # combine the data into a single data frame
  df <-
    data.frame(
      value = stan_out_hier[, i],
      group = "hier") |> 
    rbind(data.frame(
      value = stan_out_sep[, i],
      group = "sep")) |> 
    mutate(group = factor(group, levels = c("hier", "sep")))

  parts <- strsplit(i, "\\.|_")[[1]]
  
  plot_list[[i]] <- 
    ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
    # geom_density(alpha = 0.3) +
    labs(title = i, x = "Value", y = "Density") +
    geom_vline(xintercept = stan_out_hier_true[parts[2], parts[1]], linetype = "dashed", linewidth = 1, col = "red") +
    geom_vline(xintercept = stan_out_sep_true[parts[2], parts[1]], linetype = "dashed", linewidth = 1) +
    theme_minimal() +
    # cure fraction priors
    geom_density(data = cf_data, aes(x = cf_sep), alpha = 0.3, inherit.aes = FALSE, col = "black", linewidth = 1.1) +
    xlim(0, 0.75)
}

do.call(gridExtra::grid.arrange, c(plot_list, nrow = 3, ncol = data$n_endpoints))

# posterior of hierarchical cure fraction sd
sd_data <- data.frame(sd = sd)  # prior

ggplot(data = stan_out_hier, aes(x = sd_cf.1.)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  geom_density(alpha = 0.7) +
  geom_density(data = sd_data, aes(x = sd), alpha = 0.7, inherit.aes = FALSE, col = "red") +
  theme_minimal() +
  xlim(0, 2)


# plot by target statistic
plot_list <- plot_list[sort(cf_names)]

x <- do.call(gridExtra::grid.arrange, c(plot_list[1:data$n_endpoints], nrow = 3))
ggsave(plot = x, filename = glue::glue("plots/cure_fraction_posterior_ne{data$n_endpoints}.png"),
       width = 10, height = 10, bg = "white", units = "in", dpi = 300)


do.call(gridExtra::grid.arrange, c(plot_list[11:20], nrow = 3))
do.call(gridExtra::grid.arrange, c(plot_list[21:30], nrow = 3))

