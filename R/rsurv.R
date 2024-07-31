
#' Generate time to event sample with censoring
#'
#' @examples
#'
#' rsurv(n = 100,
#'       distn = "exp",
#'       prop_cens = 0.1,
#'       params = list(rate = 1))
#'
rsurv <- function(n = 100,
                  distn = "exp",
                  prop_cens = 0.1,
                  params = list(rate = 1)) {
  ##TODO: replace below  
  # rdistn <- glue("r{run$family_latent_true}")
  # latent_params_true <- eval(parse(text = run$latent_params_true))
  # rargs <- c(params, n = n)
  # times <- do.call(rdistn, rargs)
  
  t_latent <-
    if (distn == "exp") {
      do.call(rexp, c(n = n, params))
    } else if (distn == "weibull") {
      do.call(rweibull, c(n = n, params))
    } else if (distn == "biweibull") {
      do.call(rbiweibull, c(n = n, params))
    }
  
  # uniform sample then censor selected
  cens_idx <- sample(1:n,
                     size = n*prop_cens,
                     replace = FALSE)
  
  times <- t_latent
  times[cens_idx] <-
    map_dbl(times[cens_idx], function(x) runif(1, 0, x))
  
  status <- as.numeric(!(1:n) %in% cens_idx)
  
  data.frame(t_latent = t_latent,
             times = times,
             status = status)
}


#' Survival mixture cure fraction model simulation
#' for deterministic cure fraction and full
#' random sampling
#'
#' @importFrom purrr transpose
#'
rsurv_cf <- function(nsample = 20,
                     n_endpoints = 2,
                     t_cutpoint = 5,
                     mu_cf = 0.2,  # logit scale
                     sigma_cf,
                     cf_sample_method = "random",
                     cf_indiv = "fixed",
                     distn = "exp",
                     prop_cens = 0,
                     params =
                       list(
                         list(rate = 1),
                         list(rate = 1))) {
  
  if (length(params) < n_endpoints) {
    params <- rep(params, length.out = n_endpoints)
    message("recycling parameter values in end-point distributions")
  } else if (length(params) > n_endpoints) {
    stop("Number of parameter sets bigger than number of end-point distributions",
         call. = FALSE)}
  
  if (cf_sample_method == "random") {
    # hierarchically sample cure fraction
    cf_lin <- rnorm(n = n_endpoints, mu_cf, sd = sigma_cf)
    cf <- exp(cf_lin)/(1 + exp(cf_lin))
  } else if (cf_sample_method == "quantiles") {
    # generate evenly spaced probabilities
    probs <- seq(0, 1, length.out = n_endpoints + 2)[-c(1, n_endpoints + 2)]
    cf_lin <- qnorm(p = probs, mean = mu_cf, sd = sigma_cf)
    cf <- exp(cf_lin)/(1 + exp(cf_lin))
  }
    
  # # just assume the same for all endpoints to start with 
  # curestatus <- rbinom(nsample, size = 1, prob = mu_cf) + 1  # cure group indicator
  
  res <- list()
  median_times <- vector(length = n_endpoints)
  rmst <- vector(length = n_endpoints)
  
  # sample times for each endpoint
  for (i in seq_len(n_endpoints)) {
    
    res[[i]] <-
      rsurv(n = nsample,
            params = params[[i]],
            distn = distn,
            prop_cens = prop_cens)
    
    # summary statistics target measures
    median_fn <- glue("{distn}_median_cf")
    median_times[i] <- do.call(median_fn, c(params[[i]], cf = cf[i]))
    
    rmst_fn <- glue("{distn}_rmst_cf")
    rmst[i] <- do.call(rmst_fn, c(params[[i]], tmax = t_cutpoint, cf = cf[i]))
    
    if (cf_indiv == "random") {
      # sample cure status for each individual
      curestatus <- rbinom(nsample, size = 1, prob = cf[i]) + 1
    } else if (cf_indiv == "fixed") {
      # deterministic fixed size of cured
      num_cf <- round(nsample * cf[i])
      cured_idx <- sample(1:nsample, size = num_cf, replace = FALSE)
      curestatus <- rep(1, nsample)
      curestatus[cured_idx] <- 2
    }
      
    # modify times
    res[[i]] <- res[[i]] |> 
      mutate(curestatus = curestatus,
             cure_obs = curestatus == 2 & status == 1,
             after_cut = times > t_cutpoint,
             # cured
             t_latent = ifelse(cure_obs,
                               yes = NA, no = t_latent),
             times = ifelse(cure_obs,
                            yes = t_cutpoint, no = times),
             status = ifelse(cure_obs,
                             yes = 0, no = status),
             # after cut-point
             status = ifelse(after_cut,
                             yes = 0, no = status),
             times = ifelse(after_cut,
                            yes = t_cutpoint, no = times),
             endpoint = i) |> 
      select(-cure_obs, -after_cut)
  }
  
  structure(
    do.call(rbind, res),
    rmst = rmst,
    median = median_times,
    cf = cf)
}

