# simulation study for manuscript


library(multimcm)
library(dplyr)
library(survival)
library(purrr)
library(glue)

# read in scenario data
scenario_data <- read.csv(here::here("raw-data/scenarios.csv")) |> as_tibble()

# data generation
run <- scenario_data[1,]

n <- run$nsample


# sample event times
rdistn <- glue("r{run$family_latent_true}")
latent_params_true <- eval(parse(text = run$latent_params_true))

times <- NULL
for (i in 1:run$n_endpoints) {
  rargs <- c(latent_params_true[[i]], n = n)
  times <- rbind(times,
                 data.frame(endpoint = i,
                            time = do.call(rdistn, rargs)))
}


dat <- data.frame(times = rexp(n = 1000,1),
                  status = 1,
                  group = 1)

fit <- survfit(Surv(times, status) ~ group, data = dat)
plot(fit)


##TODO
# from previous analysis
# use this?
  
  
input_data <-
  rsurv_mix(n = nsample,
            n_endpoints = n_endpoints,
            t_cutpoint = t_cutpoint,
            mu_cf = mu_cf,
            sigma_cf = sigma_true,
            distn = family_latent_true,
            prop_cens = prop_censoring,
            params = latent_params_true)


########
# plots
# check Kaplan-Meier

fit_pfs <- survfit(Surv(fake_nivo_pfs$t_cens, fake_nivo_pfs$status) ~ fake_nivo_pfs$group)
fit_mix <- survfit(Surv(fake_nivo_pfs$t_cens, fake_nivo_pfs$status) ~  1)

plot(fit_pfs, xlim = c(0, 60))
lines(fit_mix, col = "blue")

  
  
  
  
  
    
############
# fit model

out <-
  bmcm_stan(
    input_data = input_data,
    formula = "Surv(time=time, event=status) ~ 1",
    cureformula = "~ treat + (1 | endpoint)",
    family_latent = family_latent,
    centre_coefs = TRUE,
    bg_model = "bg_fixed",
    bg_varname = "rate",
    bg_hr = 1,
    t_max = 400)

