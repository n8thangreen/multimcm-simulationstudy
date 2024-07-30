# performance measures for
# summary estimates for full probabilistic sampling

# see
# Using simulation studies to evaluate statistical methods (2017) Tim Morris et al, Stats in Medicine

# targets/estimates:
# * RMST of separate curves
# * cure fractions
# * median survival times

# performance measures:
# * bias
# * empirical SE
# * coverage
# * MSE

library(rstan)
library(purrr)
library(dplyr)
library(ggplot2)


scenario_data <-
  read.csv(here::here("raw-data/scenarios.csv")) |>
  as_tibble() 

data <- scenario_data[1, ]
file_name <- glue::glue("N{data$nsample}_ne{data$n_endpoints}_pcens{data$prop_censoring}_sigma{data$sigma_true}")
# data scenarios
load(paste0("data/", file_name, ".RData"))
# stan output
load(paste0("data/stan_out_", file_name, ".RData"))

load("data/input_data.RData")

target_names <- c("rmst", "median", "cf")
pm <- list()

# scenarios
# for (i in seq_along(stan_out)) {
  i <- 2  
  
  # for (j in target_names) {
    j <- "rmst"
    true_vals <- attr(input_data[[i]], which = j)
    
    pm[[j]] <- bmcm_performance_measures_N(stan_out, par_nm = j, true_vals)
  # }
# }

pm

# save(pm, file = glue::glue("data/performance_measures{file_append}.RData"))


########
# plots

load(glue::glue("data/performance_measures{file_append}.RData"))

target <- "rmst"
# target <- "cf"
# target <- "median"

# quo_measure <- quo(empirical_se)
quo_measure <- quo(bias)

plot_dat <- pm |> 
  map(target) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario") |> 
  mutate(endpoint = factor(endpoint))
  # arrange(val)

# lollipop plot

label_data <- scenario_data |> 
  mutate(label = glue::glue("{data_id}: e={n_endpoints} n={nsample} p={prop_censoring} s={sigma_true}")) |> 
  rename(scenario = id) |> 
  mutate(label = factor(label, levels = unique(label))) |> 
  select(label, scenario)

plot_dat <- plot_dat |>
  left_join(label_data, by = "scenario")

# clean measure string
measure <- quo_name(quo_measure)
y_label <- stringr::str_to_sentence(gsub(pattern = "\\_", " ", measure))

plot_dat |> 
  ggplot(aes(x = endpoint, y = !!quo_measure)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point(size=4, color="black") +
  # facet_wrap(vars(label)) +
  facet_wrap(vars(scenario)) +  # without titles
  coord_flip() +
  theme_bw() +
  xlab("Endpoint ID") +
  ylab(y_label) #+
  # ylim(0, ifelse(target == "rmst", 10, 
  #                ifelse(target == "median" & measure == "empirical_se", 5, NA)))

ggsave(filename = glue::glue("plots/lollipop_{target}_{quo_name(quo_measure)}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

# all performance measures on single plot

plot_long <- plot_dat |> 
  as.data.frame() |> 
  reshape2::melt(id.vars = c("scenario","endpoint"),
                 measure.vars = c("bias","empirical_se","mse"))

plot_long |> 
  ggplot(aes(x = endpoint, y = value)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point(size=4, color="black") +
  facet_grid(rows = vars(scenario), cols = vars(variable)) +
  coord_flip() +
  theme_bw() +
  xlab("Endpoint ID")

ggsave(filename = glue::glue("plots/lollipop_cf_priors_{target}_data_scenario_1.png"),
       height = 20, width = 20, dpi = 640, units = "cm")


#########
# tables

#' take average across curves
#'
scenario_mean_table <- function(data_list) {
  
  result <- data.frame(Scenario = character(),
                       Coverage = numeric(),
                       Bias = numeric(),
                       EmpiricalSE = numeric(),
                       stringsAsFactors = FALSE)
  
  for (i in seq_along(data_list)) {
      temp_data <- data_list[[i]]
      coverage <- mean(temp_data[, "coverage"])
      abs_bias_mean <- mean(abs(temp_data[, "bias"]))
      empirical_se <- mean(temp_data[, "empirical_se"])
      
      result <- rbind(result,
                      data.frame(Scenario = i,
                                 Coverage = coverage,
                                 Bias = abs_bias_mean,
                                 EmpiricalSE = empirical_se))
  }
  
  result
}

rmst_tab <- 
  map(pm, "rmst") |> 
  scenario_mean_table()

cf_tab <- 
  map(pm, "cf") |> 
  scenario_mean_table()

median_tab <- 
  map(pm, "median") |> 
  scenario_mean_table()

combined_tab <- 
  merge(rmst_tab, cf_tab, by = "Scenario") |> 
  merge(median_tab, by = "Scenario") |> 
  mutate(across(where(is.numeric), ~ round(., 2)))

# xtable::xtable(combined_tab, digits = 3)


library(kableExtra)

input_table <-
  scenario_data |> 
  select(id, n_endpoints, nsample, prop_censoring, sigma_true) |>
  filter(id %in% 1:16) |>
  kable(format = "latex", booktabs = TRUE)

output_table <-
  kable(combined_tab, format = "latex", booktabs = TRUE) |> 
  add_header_above(c(" ", "RMST" = 3, "Cure Fraction" = 3, "Median" = 3)) 
  # kable_styling(latex_options = c("striped", "hold_position"))

# save
write.csv(rmst_tab, file = "output_data/rmst_table.csv", row.names = FALSE)
write.csv(cf_tab, file = "output_data/cf_table.csv", row.names = FALSE)
write.csv(median_tab, file = "output_data/median_table.csv", row.names = FALSE)

