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
library(dplyr)
library(ggplot2)

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

load("data/performance_measures.RData")

# target <- "rmst"
# target <- "cf"
target <- "median"

quo_measure <- quo(empirical_se)
# quo_measure <- quo(bias)

plot_dat <- pm |> 
  map(target) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario") |> 
  mutate(endpoint = factor(endpoint))
  # arrange(val)

## drop uniquely large value  
# plot_dat <- plot_dat[-91,]
# plot_dat <- plot_dat[-11,]

# lollipop plot

# clean measure string
y_label <- stringr::str_to_sentence(gsub(pattern = "\\_", " ", quo_name(quo_measure)))

plot_dat |> 
  ggplot(aes(x = endpoint, y = !!quo_measure)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point(size=4, color="black") +
  facet_wrap(vars(scenario)) +
  coord_flip() +
  theme_bw() +
  xlab("Endpoint ID") +
  ylab(y_label)
  # ylim(0,5)

ggsave(filename = glue::glue("plots/lollipop_{target}_{quo_name(quo_measure)}.png"))


#########
# tables
