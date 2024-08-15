# to run on UCL myriad cluster
#
# script to run analysis using cmdstanr

# modify the R library path so can pick up packages that are installed in own space
.libPaths(c('/lustre/home/sejjng1/R/x86_64-pc-linux-gnu-library/4.2', .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

## install packages
# install.packages("dplyr")
# install.packages("survival")
# install.packages("purrr")
# install.packages("glue")
# install.packages("ggplot2")
# install.packages("tibble")

install.packages("epicontacts")
install.packages("tidybayes")
install.packages("here")
remotes::install_github("StatisticsHealthEconomics/multimcm", force = TRUE)

install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev',
                                       'https://cran.ma.imperial.ac.uk/', # UK
                                       getOption("repos")))
# attach
library(dplyr)
library(purrr)
library(glue)
library(tibble)
library(here)
library(parallel)
library(cmdstanr)
library(multimcm)

cmdstanr::install_cmdstan()

# functions
source("/home/sejjng1/Scratch/bmcm/functions/rsurv.R")
source("/home/sejjng1/Scratch/bmcm/functions/run_scenario.R")
source("/home/sejjng1/Scratch/bmcm/functions/target_distns.R")

# read in scenario data
scenario_data <- read.csv("/home/sejjng1/Scratch/bmcm/scenarios.csv") |> as_tibble()

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
  multimcm::precompile_bmcm_model(
    input_data = dummy_sample,
    family_latent = bmcm_params$family_latent,
    cureformula = bmcm_params$cureformula,
    use_cmdstanr = TRUE,
    file_path = "/home/sejjng1/Scratch/")

bmcm_params$precompiled_model_path <- stan_model$exe_path
  # glue::glue("/home/sejjng1/Scratch/{stan_model$model_name()}.exe")

# run for single scenario

run_scenario(1, sim_params, bmcm_params)

# parallel over multiple cores

num_cores <- 2

cl <- makeCluster(num_cores, type = "SOCK", outfile = "")

# load packages on each R process
clusterEvalQ (cl, library(cmdstanr))

# lapply(1:n_sim, \(x) run_scenario(x, sim_params, bmcm_params))
output <- parLapply(cl, 1:n_sim, fitter, sim_params, bmcm_params)

stopCluster(cl)





