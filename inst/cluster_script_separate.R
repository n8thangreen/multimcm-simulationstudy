# to run on UCL myriad cluster
#
# script to run analysis using cmdstanr
# for separate model

# get task ID from the command line arguments
args <- commandArgs(trailingOnly = TRUE)

task_id <- as.integer(args[1])
SCENARIO_ID <- as.integer(args[2])

# modify the R library path so can pick up packages that are installed in own space
.libPaths(c('/lustre/home/sejjng1/R/x86_64-pc-linux-gnu-library/4.2', .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

packages <- c("epicontacts", "tidybayes", "here")
installed_packages <- rownames(installed.packages())

for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

## install packages

if (!("multimcm" %in% installed_packages)) {
  remotes::install_github("StatisticsHealthEconomics/multimcm")
}

if (!("cmdstanr" %in% installed_packages)) {
  install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev',
                                         'https://cran.ma.imperial.ac.uk/', # UK
                                         getOption("repos")))
}

# attach
library(dplyr)
library(purrr)
library(glue)
library(tibble)
library(here)
library(cmdstanr)
library(parallel)
library(multimcm)

if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  cmdstanr::install_cmdstan()
}

# functions
source("/home/sejjng1/Scratch/bmcm/functions/rsurv.R")
source("/home/sejjng1/Scratch/bmcm/functions/run_scenario.R")
source("/home/sejjng1/Scratch/bmcm/functions/target_distns.R")

# read in scenario data
scenario_data <- read.csv(here::here("/home/sejjng1/Scratch/bmcm/scenarios.csv")) |> as_tibble()

# scenario number
i <- SCENARIO_ID

data <- scenario_data[i, ]
latent_params_true <- eval(parse(text = data$latent_params_true))

# prior parameters for each cure fraction
# duplicate for each treatment
prior_cure_list <- 
  rep(list(mu_alpha = rep(data$mu_cf_prior, 2),
       sigma_alpha = rep(data$sigma_cf_prior, 2)), data$n_endpoints)

names(prior_cure_list) <-
  paste0(names(prior_cure_list), "_", rep(1:data$n_endpoints, each = 2))

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
    cureformula = "~ tx + endpoint",   # separate model
    family_latent = data$family_latent_model,
    prior_cure = prior_cure_list,   
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

bmcm_params$precompiled_model_path <- stan_model$exe_file()
# glue::glue("/home/sejjng1/Scratch/{stan_model$model_name()}.exe")

# check if scenario output folder exits and create
output_dir <- glue::glue("/home/sejjng1/Scratch/output/scenario_{data$data_id}")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# scenario run

max_attempts <- 2
attempt <- 0
success <- FALSE

while (!success && attempt < max_attempts) {
  attempt <- attempt + 1
  paste("Attempt", attempt, "of", max_attempts, "\n")
  
  fit <- try(
    run_scenario(task_id, sim_params, bmcm_params, dir = glue::glue("{output_dir}/"))
  )
  
  # check model ran successfully
  if (inherits(fit, "CmdStanMCMC")) success <- TRUE
}



