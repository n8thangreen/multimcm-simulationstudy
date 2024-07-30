# fit Bayesian mixture cure models
# using simulation data


library(parallel)
library(doSNOW)
library(dplyr)
library(multimcm)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

# load simulation data
# load(file = "data/determ_input_data.RData")

# full probabilistic for all scenarios
datasets_all <- vector(mode = "list", 16)
for (i in 1:16) {
  data <- scenario_data[i, ]
  file_name <- glue::glue("N{data$nsample}_ne{data$n_endpoints}_pcens{data$prop_censoring}_sigma{data$sigma_true}")
  load(paste0("data/", file_name, ".RData"))
  datasets_all[[i]] <- input_data
}

###############################
# full probabilistic sampling

stan_out <- list()

i <- 1
n_sim <- length(datasets_all[[1]])
input_dat <- list()
  
# for (i in 1:nrow(scenario_data)) {
    
  params <- scenario_data[i, ]
    
  # duplicate for each treatment
  prior_cure <- list(mu_alpha = rep(params$mu_cf_prior, 2),
                     sigma_alpha = rep(params$sigma_cf_prior, 2),
                     mu_sd_cf = rep(params$mu_sd_cf_prior, 2),
                     sigma_sd_cf = rep(params$sigma_sd_cf_prior, 2))
    
  stan_out[[i]] <- list()

  # include additional variables in dataset
  # need for bmcm_stan()
  
  input_dat[[i]] <- list()

  for (j in 1:n_sim) {
    input_dat[[i]][[j]] <- mutate(datasets_all[[i]][[j]],
                                  tx = 1,           # only a single treatment
                                  rate = 10^(-10))  # background hazard
    
    ##TODO; errors with single treatment only
    ## quick fix by duplicating inputs
    input_dat[[i]][[j]] <- rbind(input_dat[[i]][[j]],
                                 mutate(datasets_all[[i]][[j]],
                                        tx = 2,
                                        rate = 10^(-100)))
  }
    
  stan_out[[i]] <-
    lapply(input_dat[[i]], function(x) {
      bmcm_stan(
        input_data = x,
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
    })
# }

# save(stan_out, file = "data/stan_out.RData")
  
# using parallel package
num_cores <- detectCores() - 1

cl <- makeCluster(num_cores, type = "SOCK", outfile = "")

clusterEvalQ(cl, {
  library(ggplot2)
  library(stringr)
})

res <- parLapply(
  cl, input_dat[[i]],
  bmcm_stan,
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

stopCluster(cl)


## using foreach

cl <- makeCluster(num_cores, type = "SOCK", outfile = "")
registerDoSNOW(cl)

pb <- txtProgressBar(max = n_sim, style = 3)

progress <- function(n)
  setTxtProgressBar(pb, n)

opts <- list(progress = progress)

res <- foreach(
  j = 1:n_sim,
  # .combine = 'comb',
  # .multicombine = TRUE,
  # .init = list(list(), list(), list(), list()),
  .options.snow = opts,
  .packages = c("boot")
) %dopar% {
  bmcm_stan(input_dat[[i]][[j]],
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

stopCluster(cl)

  
##############################
# deterministic cure fraction

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

save(stan_out, file = "data/stan_out.RData")

#######
# plot

for (i in seq_along(stan_out)) {
  gg <- plot_S_joint(stan_out[[i]], add_km = TRUE) + xlim(0,5) + facet_wrap(vars(endpoint))
  ggsave(plot = gg, device = "png", filename = glue::glue("plots/survival_plots_{i}.png"))
}




