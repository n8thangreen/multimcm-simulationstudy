# tests


## large sample size and number of endpoint

input_data <-
  rsurv_mix(n = 1000,
            n_endpoints = 20,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "weibull",
            prop_cens = 0,
            params = list(list(shape = 1, scale = 1)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

## large sample size and small number of endpoint

input_data <-
  rsurv_mix(n = 1000,
            n_endpoints = 3,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "weibull",
            prop_cens = 0,
            params = list(list(shape = 1, scale = 1)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

## small sample size and small number of endpoint

input_data <-
  rsurv_mix(n = 20,
            n_endpoints = 3,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "weibull",
            prop_cens = 0,
            params = list(list(shape = 1, scale = 1)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

## small sample size and large number of endpoint

input_data <-
  rsurv_mix(n = 20,
            n_endpoints = 20,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "weibull",
            prop_cens = 0,
            params = list(list(shape = 1, scale = 1)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")


## different rates

input_data <-
  rsurv_mix(n = 1000,
            n_endpoints = 20,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "exp",
            prop_cens = 0,
            params = list(list(rate = 10)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")


# two distributions
input_data <-
  rsurv_mix(n = 1000,
            n_endpoints = 20,
            t_cutpoint = 6,
            mu_cf = -1.39,
            sigma_cf = 2,
            cf_sample_method = "quantiles",
            distn = "exp",
            prop_cens = 0,
            params = list(list(rate = 1),
                          list(rate = 0.1)))

# check Kaplan-Meier

# single plots
fit_dat <- survfit(Surv(times, status) ~ endpoint, data = input_data)

plot(fit_dat, col = rgb(0,0,0, alpha = 0.1))
abline(h = exp(-1.39)/(1 + exp(-1.39)), col = "red")

