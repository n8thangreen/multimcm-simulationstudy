# simulation study for manuscript


library(multimcm)
library(dplyr)
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


  
  

##TODO
# from previous analysis
# use this?
  
  
s_tx_groups <- table(surv_input_data$TRTA)
mu0 <- -2
OSage_NIVO <- surv_input_data$OSage[surv_input_data$TRTA == "NIVOLUMAB"]

fake_nivo_pfs <-
  rsurv_mix(cf = 0.45,
            n = s_tx_groups["NIVOLUMAB"],
            distn = c("exp", "exp"),
            prop_cens = 0.2,
            params =
              list(
                list(
                  mu = c(mu0, 0.005)),
                list(
                  mu = c(-8.5, 0.05))),
            X = OSage_NIVO)

# centered for joint model
surv_input_fake <-
  surv_input_data %>%
  mutate(pfs = ifelse(TRTA == "NIVOLUMAB",
                      fake_nivo_pfs$t_cens, pfs),
         pfs_event = ifelse(TRTA == "NIVOLUMAB",
                            fake_nivo_pfs$status, pfs_event)) %>%
  group_by(TRTA) %>%
  mutate(pfs_centred = pfs - 1/exp(mu0))

X_nivo <-
  data.frame(
    OSage = OSage_NIVO,
    t_pfs = surv_input_fake$pfs_centred[surv_input_fake$TRTA == "NIVOLUMAB"])

fake_nivo_os <-
  rsurv_mix(cf = 0.2,
            n = s_tx_groups["NIVOLUMAB"],
            distn = c("exp", "exp"),
            prop_cens = 0.2,
            params =
              list(
                list(
                  mu = c(-3, 0.005, -0.001)),
                list(
                  mu = c(-8.5, 0.03, 0))),
            X = X_nivo)

# replace with fake data
surv_input_fake$os[surv_input_fake$TRTA == "NIVOLUMAB"] <- fake_nivo_os$t_cens
surv_input_fake$os_event[surv_input_fake$TRTA == "NIVOLUMAB"] <- fake_nivo_os$status


## plots
# check Kaplan-Meier

library(survival)
fit_pfs <- survfit(Surv(fake_nivo_pfs$t_cens, fake_nivo_pfs$status) ~ fake_nivo_pfs$group)
fit_mix <- survfit(Surv(fake_nivo_pfs$t_cens, fake_nivo_pfs$status) ~  1)
plot(fit_pfs, xlim = c(0, 60))
lines(fit_mix, col = "blue")

fit_os <- survfit(Surv(fake_nivo_os$t_cens, fake_nivo_os$status) ~ fake_nivo_os$group)
fit_mix <- survfit(Surv(fake_nivo_os$t_cens, fake_nivo_os$status) ~  1)
plot(fit_os, xlim = c(0, 60))
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

