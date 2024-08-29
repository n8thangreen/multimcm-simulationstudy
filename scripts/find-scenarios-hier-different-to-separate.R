# find scenarios where the separate model is different to the hierarchical model
# then use this in the simulation study


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)
library(ggplot2)


###########
# simulate 

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

# cure fraction
## true
invlogit(rnorm(1000, -1.39, 0.1)) |> 
  density() |> plot()

## prior
n <- 10000
mu <- rnorm(n, mean = -1.39, sd = 0.1)
mu_cf <- invlogit(mu)
sd <- truncated_cauchy(n, location = 0.15, scale = 0.1, 0)
sd <- sd[sd<10]

density(mu) |> plot()
density(mu_cf) |> plot()
density(sd) |> plot()

invlogit(rnorm(n = n, mean = mu, sd = sd)) |> 
  density() |> plot(xlim = c(0,0.5))


# latent survival
n <- 100
shape <- rgamma(n, shape = 2, scale = 1)
scale <- rlnorm(n, meanlog = 0.6, sdlog = 0.5)

density(shape) |> plot()
density(scale) |> plot()

plot(NULL, ylim = c(0,1), xlim = c(0,5))
for (i in 1:n) {
  pweibull(q = seq(0,5,0.1), shape = shape[i], scale = scale[i], lower.tail = F) |> 
    lines(x = seq(0,5,0.1), col = "grey")
}

# true survival curves
pweibull(q = seq(0,5,0.1), shape = 1, scale = 4, lower.tail = F) |> 
    plot(x = seq(0,5,0.1), col = "red", xlim = c(0,5), ylim = c(0,1), type = "l")
pweibull(q = seq(0,5,0.1), shape = 1, scale = 1, lower.tail = F) |> 
    lines(x = seq(0,5,0.1))

#############
# analysis

# scenario_data
data <- data.frame(
  nsample = 10,
  n_endpoints = 10,
  t_cutpoint = 5,
  mu_cf_prior = -1.39,
  sigma_cf_prior = 0.1,
  mu_sd_cf_prior = 0.15,
  sigma_sd_cf_prior = 0.1,
  cf_true = -1.39,
  sigma_true = 0.1,
  family_latent_true = "weibull",
  family_latent_model = "weibull",
  prop_censoring = 0,
  latent_params_true =
    "list(list(shape = 1, scale = 1))",
    # "list(list(shape = 1, scale = 4), list(shape = 1, scale = 1))",  # scale is 1/rate
  a_shape_latent_prior = 2,       # gamma on shape
  b_shape_latent_prior = 1,
  mu_S_prior = 0.6,               # log-normal on scale
  sigma_S_prior = 0.5)  

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

# prior parameters for each latent curve
# duplicate for each end point
prior_latent_list <- 
  rep(list(a_shape = data$a_shape_latent_prior,
           b_shape = data$b_shape_latent_prior,
           mu_S = as.array(data$mu_S_prior),
           sigma_S = as.array(data$sigma_S_prior)), data$n_endpoints)
names(prior_latent_list) <-
  paste0(names(prior_latent_list), "_", rep(1:data$n_endpoints, each = 4))

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

run_scenario(1, sim_params, bmcm_params_hier, dir = "output_data/hierarchical/", rstan_format = FALSE, seed = 1234)
run_scenario(1, sim_params, bmcm_params_sep, dir = "output_data/separate/", rstan_format = FALSE, seed = 1234)


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

plot_list <- list()

for (i in cf_names) {
  # combine the data into a single data frame
  df <- data.frame(
    value = c(stan_out_hier[, i], stan_out_sep[, i]),
    group = factor(rep(c("hier", "sep"), each = 500))
  )
  
  parts <- strsplit(i, "\\.|_")[[1]]
  
  plot_list[[i]] <- 
    ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
    geom_density(alpha = 0.7) +
    labs(title = i, x = "Value", y = "Density") +
    geom_vline(xintercept = stan_out_hier_true[parts[2], parts[1]], linetype = "dashed", linewidth = 1, col = "red") +
    geom_vline(xintercept = stan_out_sep_true[parts[2], parts[1]], linetype = "dashed", linewidth = 1) +
    theme_minimal()
}

do.call(gridExtra::grid.arrange, c(plot_list, nrow = 3, ncol = data$n_endpoints))

# separate measures
do.call(gridExtra::grid.arrange, c(plot_list[1:data$n_endpoints], nrow = 3, ncol = ceiling(data$n_endpoints/3)))
do.call(gridExtra::grid.arrange, c(plot_list[(data$n_endpoints+1):(2*data$n_endpoints)], nrow = 3, ncol = ceiling(data$n_endpoints/3)))
do.call(gridExtra::grid.arrange, c(plot_list[(2*data$n_endpoints+1):(3*data$n_endpoints)], nrow = 3, ncol = ceiling(data$n_endpoints/3)))
