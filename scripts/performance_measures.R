# performance measures for
# summary estimates

# see
# Using simulation studies to evaluate statistical methods (2017) Tim Morris et al, Stats in Medicine

##TODO: 
## * take average across curves?

# targets/estimates:
# * RMST of separate curves
# * cure fractions
# * median survival times

# performance measures:
# * bias
# * empirical SE
# * coverage
# * MSE


##TODO: how to use posterior package?

library(rstan)
library(purrr)

load("data/stan_out.RData")
load("data/input_data.RData")

target_names <- c("rmst", "median", "cf")
pm <- list()

# scenarios
for (i in seq_along(stan_out)) {
  fit <- stan_out[[i]]
  pm[[i]] <- list()
  
  # estimates
  for (j in target_names) {
    true_vals <- attr(input_data[[i]], which = j)
    
    pm[[i]][[j]] <- bmcm_performance_measures(fit, par_nm = j, true_vals)
  }
}

pm

save(pm, file = "data/performance_measures.RData")


########
# plots

# measure <- "rmst"
# measure <- "cf"
measure <- "median"

plot_dat <- pm |> 
  map(measure) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario")

## drop uniquely large value  
plot_dat <- plot_dat[-91,]

# lollipop plot

# bias

plot_dat %>%
  # arrange(val) |> 
  # mutate(name=factor(name, levels=name)) |>    # update the factor levels
  ggplot(aes(x = endpoint, y = bias)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point( size=4, color="black") +
  facet_wrap(vars(scenario)) +
  coord_flip() +
  theme_bw() +
  xlab("")

# se
## drop uniquely large value  
plot_dat <- plot_dat[-11,]

plot_dat %>%
  ggplot(aes(x = endpoint, y = empirical_se)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point( size=4, color="black") +
  facet_wrap(vars(scenario)) +
  coord_flip() +
  theme_bw() +
  xlab("")


#########
# tables
