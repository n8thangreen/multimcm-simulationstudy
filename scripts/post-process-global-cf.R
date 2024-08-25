# performance measures for
# summary estimates
# global cure fraction

library(rstan)
library(purrr)
library(dplyr)
library(ggplot2)


scenario_data <-
  read.csv(here::here("raw-data/scenarios.csv")) |>
  as_tibble() 

# # data scenarios
# file_append <- ""

# prior scenarios
file_append <- "_cf_priors"

load(glue::glue("data/stan_out{file_append}.RData"))

load("data/input_data.RData")

#######################
# performance measures

true_val_names <- c("cf_true", "sigma_true")

# replicate the first scenario values
true_vals <- scenario_data[1, true_val_names] |> 
  rename(lp_cf_global = cf_true, sd_cf = sigma_true) |> 
  # slice(rep(1:n(), each = 5)) |> 
  as.list()

pm_cf <- list()

for (i in seq_along(stan_out)) {
  fit <- stan_out[[i]]
  pm_cf[[i]] <- list()
  
  # estimates
  for (j in c("lp_cf_global", "sd_cf")) {
    
    pm_cf[[i]][[j]] <- bmcm_performance_measures(fit, par_nm = j, true_vals[[j]], append = FALSE)
  }
}

########
# plots

library(tidyr)
library(bayesplot)

# posteriors

# draw grid of plot using ggplot2 with true value marked by vertical line
# for each cure fraction scenario in stan_out

plot_dat <- map(stan_out, "output")

stan_extract <- map(plot_dat, rstan::extract, pars = "sd_cf")
xx <- map(stan_extract, ~as.matrix(.x[[1]])[,1])
xx <- do.call(cbind, xx)
colnames(xx) <- paste0(1:ncol(xx))

p1 <- mcmc_intervals(xx, 
                     prob = 0.8,           # 80% credible intervals
                     prob_outer = 0.95) +  # 95% credible intervals
  geom_vline(xintercept = true_vals$sd_cf, linetype = "dashed", color = "red") +
  xlab("Standard deviation of global cure fraction") +
  theme_minimal()

stan_extract <- map(plot_dat, rstan::extract, pars = "lp_cf_global")
xx <- map(stan_extract, ~as.matrix(.x[[1]])[,1])
xx <- do.call(cbind, xx)
colnames(xx) <- paste0(1:ncol(xx))
p2 <- mcmc_intervals(xx, 
                     prob = 0.8,           # 80% credible intervals
                     prob_outer = 0.95) +  # 95% credible intervals
  geom_vline(xintercept = true_vals$lp_cf_global, linetype = "dashed", color = "red") +
  xlab("Mean of global cure fraction") +
  theme_minimal()

grid_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(grid_plot, filename = glue::glue("plots/posterior_forest_plot_global_cf.png"),
       height = 20, width = 20, dpi = 640, units = "cm")


load(file = "data/cure_fraction_hyperparameters.RData")


#################
# include priors

# Extract data from the existing plot
plot_data <- mcmc_intervals_data(xx, 
                                 prob = 0.8, 
                                 prob_outer = 0.95)
## mean

# additional data for new bars
mean_data <- plot_data
for (i in 1:5) {
  mean_data$ll[i] <- params$mean[[i]][[1]] - 1.96*params$mean[[i]][[2]]
  mean_data$l[i] <- params$mean[[i]][[1]] - 1.28*params$mean[[i]][[2]]
  mean_data$m[i] <- params$mean[[i]][[1]]
  mean_data$h[i] <- params$mean[[i]][[1]] + 1.28*params$mean[[i]][[2]]
  mean_data$hh[i] <- params$mean[[i]][[1]] + 1.96*params$mean[[i]][[2]]
}

height <- 0

# Add the additional bars to the plot
p3 <- p2 + 
  xlim(-4,2) +
  geom_linerange(aes(y = parameter, xmin = ll, xmax = hh), 
                 data = mean_data, 
                 height = height,    # Adjust the height to position below the main bars
                 color = "lightgreen",  # Example color for the additional bars
                 position = position_nudge(y = -0.2)) +
  geom_linerange(aes(y = parameter, xmin = l, xmax = h), 
                 data = mean_data,
                 linewidth = 2,
                 height = height,    # Adjust the height to position below the main bars
                 color = "darkgreen",  # Example color for the additional bars
                 position = position_nudge(y = -0.2)) +
  geom_point(aes(y = parameter, x = m), 
             data = mean_data, 
             shape = 21, 
             fill = "lightgreen", 
             color = "black",
             position = position_nudge(y = -0.2),
             size = 4) 
  # geom_text(aes(y = parameter, x = 1, label = glue::glue("{m} [{round(l,2)},{round(h,2)}]")), 
  #           data = mean_data, 
  #           hjust = -0.1,   # Adjust horizontal justification to position text outside the plot area
  #           size = 3.5,     # Adjust text size as needed
  #           color = "black",
  #           vjust = 0.5)
  

## sd

# additional data for new bars
sd_data <- plot_data
for (i in 1:5) {
  sd_data$ll[i] <- params$sd[[i]][[1]] - 1.96*params$sd[[i]][[2]]
  sd_data$l[i] <- params$sd[[i]][[1]] - 1.28*params$sd[[i]][[2]]
  sd_data$m[i] <- params$sd[[i]][[1]]
  sd_data$h[i] <- params$sd[[i]][[1]] + 1.28*params$sd[[i]][[2]]
  sd_data$hh[i] <- params$sd[[i]][[1]] + 1.96*params$sd[[i]][[2]]
}

library(ggrepel)

# Add the additional bars to the plot
p4 <- p1 + 
  xlim(-5,6) +
  geom_linerange(aes(y = parameter, xmin = ll, xmax = hh), 
                 data = sd_data, 
                 height = height,    # Adjust the height to position below the main bars
                 color = "lightgreen",  # Example color for the additional bars
                 position = position_nudge(y = -0.2)) +
  geom_linerange(aes(y = parameter, xmin = l, xmax = h), 
                 data = sd_data,
                 linewidth = 2,
                 height = height,    # Adjust the height to position below the main bars
                 color = "darkgreen",  # Example color for the additional bars
                 position = position_nudge(y = -0.2)) +
  geom_point(aes(y = parameter, x = m), 
             data = sd_data, 
             shape = 21, 
             fill = "lightgreen", 
             color = "black",
             position = position_nudge(y = -0.2),
             size = 4) 
  # geom_text(aes(y = parameter, x = 2.5, label = glue::glue("{m} [{round(l,2)},{round(h,2)}]")), 
  #           data = sd_data, 
  #           hjust = -0.1,   # Adjust horizontal justification to position text outside the plot area
  #           size = 3.5,     # Adjust text size as needed
  #           color = "black",
  #           vjust = 0.5) 

grid_plot <- gridExtra::grid.arrange(p3, p4, ncol = 2)

ggsave(grid_plot, filename = glue::glue("plots/posterior_forest_plot_global_cf.png"),
       height = 20, width = 20, dpi = 640, units = "cm")
