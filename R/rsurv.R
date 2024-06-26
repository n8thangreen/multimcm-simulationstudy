
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
  times <-
    if (distn == "exp") {
      rexp(n, params$rate)
    } else if (distn == "weibull") {
      rweibull(n, params$alpha, params$mu)
    } else if (distn == "biweibull") {
      do.call(rbiweibull, c(n = n, params))
    }
  
  # uniform sample then censor selected
  cens_idx <- sample(1:n,
                     size = n*prop_cens,
                     replace = FALSE)
  
  t_cens <- times
  t_cens[cens_idx] <-
    map_dbl(t_cens[cens_idx], function(x) runif(1, 0, x))
  
  status <- as.numeric(!(1:n) %in% cens_idx)
  
  data.frame(times = times,
             t_cens = t_cens,
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
  cf <- rnorm(n = n_endpoints, mu_cf, sd = sigma_cf)
  
  # just assume the same for all endpoints to start with 
  curestatus <- rbinom(nsample, size = 1, prob = mu_cf) + 1  # cure group indicator
  ncf <- c(sum(curestatus == 1), sum(curestatus == 2))       # group sizes
  
  res <- list()
  
  # sample times for each endpoint
  for (i in seq_len(n_endpoints)) {
    
    res[[i]] <-
      rsurv(n = nsample,
            params = params[[i]],
            distn = distn,
            prop_cens = prop_cens)
    
    # # hierarchically sample cure status
    # curestatus <- rbinom(nsample, size = 1, prob = cf) + 1  # cure group indicator
    
    # modify times
    res[[i]] <- res[[i]] |> 
      mutate(curestatus = curestatus,
             # cured
             times = ifelse(curestatus == 2 & status == 1,
                            yes = t_cutpoint, no = times),
             t_cens = ifelse(curestatus == 2 & status == 1,
                             yes = t_cutpoint, no = t_cens),
             # after cut-point
             status = ifelse(t_cens > t_cutpoint,
                             yes = 0, no = status),
             t_cens = ifelse(t_cens > t_cutpoint,
                             yes = t_cutpoint, no = t_cens)
      )
  }
  
  ##TODO
  # what format to return?
  out <- purrr::transpose(res)
  
  list(
    times = unlist(out$times),
    t_cens = unlist(out$t_cens),
    status = unlist(out$status),
    curestatus = unlist(out$curestatus))
}

