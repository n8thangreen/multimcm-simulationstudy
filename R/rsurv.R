
#' Generate time to event sample with censoring
#'
#' @examples
#'
# rsurv(n = 100,
#       distn = "exp",
#       prop_cens = 0.1,
#       params = list(rate = 1))
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
#'
#' @importFrom purrr transpose
#'
rsurv_mix <- function(nsample = 20,
                      n_endpoints = 2,
                      t_cutpoint,
                      mu_cf = 0.2,
                      sigma_cf,
                      distn = "exp",
                      prop_cens = 0,
                      params =
                        list(
                          list(rate = 1),
                          list(rate = 1))) {
  
  if (length(params) != n_endpoints)
    stop("Number of parameter sets and distributions don't match",
         call. = FALSE)
  
  # hierarchically sample cure fraction
  cf_lin <- rnorm(n = n_endpoints, mu_cf, sd = sigma_cf)
  cf <- exp(cf_lin)/(1 + exp(cf_lin))
    
  # # just assume the same for all endpoints to start with 
  # curestatus <- rbinom(nsample, size = 1, prob = mu_cf) + 1  # cure group indicator
  
  res <- list()
  
  # sample times for each endpoint
  for (i in seq_len(n_endpoints)) {
    
    res[[i]] <-
      rsurv(n = nsample,
            params = params[[i]],
            distn = distn,
            prop_cens = prop_cens)
    
    # hierarchically sample cure status
    curestatus <- rbinom(nsample, size = 1, prob = cf[i]) + 1  # cure group indicator
    
    # modify times
    res[[i]] <- res[[i]] |> 
      mutate(curestatus = curestatus,
             # cured
             t_latent = ifelse(curestatus == 2 & status == 1,
                            yes = t_cutpoint, no = t_latent),
             times = ifelse(curestatus == 2 & status == 1,
                             yes = t_cutpoint, no = times),
             # after cut-point
             status = ifelse(times > t_cutpoint,
                             yes = 0, no = status),
             times = ifelse(times > t_cutpoint,
                             yes = t_cutpoint, no = times),
             endpoint = i)
  }
  
  do.call(rbind, res)
}

