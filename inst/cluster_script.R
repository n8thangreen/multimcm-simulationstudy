# to run on UCL myriad cluster
#
# script to run analysis using cmdstanr

# Get the task ID from the command line arguments
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
# install.packages("dplyr")
# install.packages("survival")
# install.packages("purrr")
# install.packages("glue")
# install.packages("ggplot2")
# install.packages("tibble")

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
i <- 2

# from submitted job with Environment Variables
exists("SCENARIO_ID")
if (exists("SCENARIO_ID")) i <- SCENARIO_ID

data <- scenario_data[i, ]
latent_params_true <- eval(parse(text = data$latent_params_true))

# number of parallel runs
#n_sim <- 2

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

bmcm_params$precompiled_model_path <- stan_model$exe_file()
# glue::glue("/home/sejjng1/Scratch/{stan_model$model_name()}.exe")

# check if scenario output folder exits and create
output_dir <- glue::glue("/home/sejjng1/Scratch/output/scenario_{data$data_id}")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# run for single scenario

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

# parallel over multiple cores
#
#num_cores <- 20
#
#cl <- makeCluster(num_cores, type = "SOCK", outfile = "")
#
# load packages on each R process
#clusterEvalQ (cl, {
#  library(cmdstanr)
#  library(purrr)
#  library(dplyr)
#  library(glue)
#  library(multimcm)
#})
#
# export functions to the cluster
#clusterExport(cl, varlist = c("rsurv_cf", "rsurv",    # functions
#			      "run_scenario",
#                              "extract_params",
#                              "weibull_rmst_cf", "weibull_rmst",
#                              "weibull_median_cf", "weibull_median",
#			      "sim_params", "bmcm_params"))  # variables
#
# lapply(1:n_sim, \(x) run_scenario(x, sim_params, bmcm_params, "/home/sejjng1/Scratch/"))
#output <- parLapply(cl, 1:n_sim, function(i) run_scenario(i, sim_params, bmcm_params, "/home/sejjng1/Scratch/"))
#
#stopCluster(cl)





