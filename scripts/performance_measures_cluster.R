# performance measures for
# summary estimates for full probabilistic sampling
# with samples data produced using the cluster
# for either separate or hierarchical models
#
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


# input values
scenario_data <-
  read.csv(here::here("raw-data/scenarios.csv")) |>
  as_tibble() 

n_scenario <- nrow(scenario_data)

target_names <- c("rmst", "median", "cf")
pm <- list()
summary_dat <- list()

## select
model_type <- "separate"
# model_type <- "hierarchical"

# for each scenario get sample summary statistics
# and performance measures
for (i in 1:n_scenario) {
  print(i)
  
  data <- scenario_data[i, ]
  
  cluster_data_dir <-
    glue::glue("output_data/cluster/mu 0.01/{model_type}/scenario_{i}/")
  
  # input data
  true_vals_filenames <- dir(path = cluster_data_dir,
                             pattern = "^true_values",
                             full.names = TRUE)
  true_vals <- lapply(true_vals_filenames, \(x) read.csv(x))
  
  # stan output
  stan_filenames <- dir(path = cluster_data_dir,
                        pattern = "^samples", full.names = TRUE)
  stan_out <- lapply(stan_filenames, \(x) read.csv(x))
  
  # clean column names
  # because treatments are the same so can combine
  stan_out <- 
    lapply(stan_out, \(x) {
      # remove the .number. from the column names
      colnames(x) <- sub("\\.\\d+\\.$", "", colnames(x))
      
      # combine same column names  
      df_long <- x |> 
        tidyr::pivot_longer(cols = -X,
                            names_to = "variable",
                            values_to = "value") |> 
        arrange(variable) |> 
        group_by(variable) |>
        mutate(id = 1:n()) |> 
        reshape2::dcast(id ~ variable, value.var = "value")
    })
  
  nendpoints <- nrow(true_vals[[1]])
  
  pm[[i]] <- list()
  summary_dat[[i]] <- list()
  
  for (j in target_names) {
    
    summary_dat[[i]][[j]] <- 
      lapply(1:nendpoints,
             \(x) samples_summary_stats(
               stan_out,
               par_nm = j,
               true_vals,
               endpoint_id = x))
    
    pm[[i]][[j]] <-
      performance_measures_cluster(stan_out, par_nm = j, true_vals)
  }
}

pm

save(pm, file = glue::glue("data/performance_measures_cluster_{model_type}.RData"))
save(summary_dat, file = glue::glue("data/summary_data_cluster_{model_type}.RData"))


########
# plots
########

load(glue::glue("data/performance_measures_cluster_{model_type}.RData"))
load(glue::glue("data/summary_data_cluster_{model_type}.RData"))

##########################
# histograms of theta_hat

target <- "rmst"
# target <- "cf"
# target <- "median"

# endp <- 1
endp <- 2

hist_dat <- 
  summary_dat |> 
  map(target) |> 
  map(~ .x[[endp]])

x_min <- min(map_dbl(hist_dat, ~ min(.x$theta_true)))
x_max <- max(map_dbl(hist_dat, ~ max(.x$theta_true)))

plot_list <- map(
  seq_along(hist_dat), 
  ~ ggplot(hist_dat[[.x]], aes(x = theta_hat)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = mean(hist_dat[[.x]]$theta_true), color = "red") +
    ggtitle(glue::glue("Scenario {.x}")) +
    xlab(target) +
    ylab("Frequency") +
    xlim(x_min, x_max) +
    theme_bw())

patchwork::wrap_plots(plot_list, ncol = 4)

ggsave(filename = glue::glue("plots/theta_hat_hist_{target}_endpoint{endp}_{model_type}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

#######################
# zip plot of coverage

# target <- "rmst"
target <- "cf"
# target <- "median"

endp <- 1
# endp <- 2

zip_dat <- 
  summary_dat |> 
  map(target) |> 
  map(~ .x[[endp]])

xx <- map(zip_dat,
          ~mutate(.x,
             z = (theta_true - theta_hat)*4/(theta_hat_upp - theta_hat_low),
             z_mean = (mean(theta_true) - theta_hat)*4/(theta_hat_upp - theta_hat_low),
             p_value = 2*pnorm(-abs(z)),
             p_value_mean = 2*pnorm(-abs(z_mean)),
             sig = p_value < 0.05,
             sig_mean = p_value_mean < 0.05)|> 
  arrange(desc(p_value)) |> 
  mutate(order = 1:n()) |> 
  arrange(desc(p_value_mean)) |> 
  mutate(order_mean = 1:n()))

# zip_pop
plot_list <- 
  map(xx,
      ~ggplot(.x) +
        geom_segment(aes(y = order, x = theta_hat_low, xend = theta_hat_upp, col = sig)) + 
        geom_point(aes(y=order, x=theta_true), color = "grey", size = 0.3, inherit.aes = F) +
        geom_hline(yintercept = 950, linewidth = 1.2, linetype = "dashed") +
        theme_bw() +
        xlab(target) +
        theme(legend.position="none"))

patchwork::wrap_plots(plot_list, nrow = 4)

ggsave(filename = glue::glue("plots/zip_pop_{target}_endpoint{endp}_{model_type}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

# zip_pop_mean
plot_list <- 
  map(xx,
      ~ggplot(.x) +
        geom_segment(aes(y = order_mean, x = theta_hat_low, xend = theta_hat_upp, col = sig_mean)) + 
        geom_vline(xintercept = mean(.x$theta_true), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 950, size = 1.2, linetype = "dashed") +
        theme_bw() +
        xlab(target) +
        theme(legend.position="none"))

patchwork::wrap_plots(plot_list, ncol = 4)

ggsave(filename = glue::glue("plots/zip_pop_mean_{target}_endpoint{endp}_{model_type}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

#################
# lollipop plots

# target <- "rmst"
target <- "cf"
# target <- "median"

# quo_measure <- quo(bias)
# quo_measure <- quo(relative_bias)
# quo_measure <- quo(coverage)
quo_measure <- quo(empirical_se)

# extract target data and combine scenarios
plot_dat <- pm |> 
  map(target) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario") |> 
  mutate(endpoint = factor(endpoint))
# arrange(val)

label_data <-
  scenario_data |> 
  mutate(label = glue::glue("{data_id}: e={n_endpoints} n={nsample}")) |> 
  # mutate(label = glue::glue("{data_id}: e={n_endpoints} n={nsample} p={prop_censoring} s={sigma_true}")) |> 
  rename(scenario = data_id) |> 
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
  facet_wrap(vars(label), nrow = 4) +
  # facet_wrap(vars(scenario)) +  # without titles
  coord_flip() +
  theme_bw() +
  xlab("Endpoint ID") +
  ylab(y_label) #+
# ylim(0, 0.5)
# ylim(0, ifelse(target == "rmst", 10, 
#                ifelse(target == "median" & measure == "empirical_se", 5, NA)))

ggsave(filename = glue::glue("plots/lollipop_{target}_{quo_name(quo_measure)}_{model_type}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

# all performance measures on single plot
# used for cure fraction prior base case scenario

plot_long <- plot_dat |> 
  as.data.frame() |> 
  reshape2::melt(id.vars = c("scenario","endpoint"),
                 measure.vars = c("bias","empirical_se","mse","coverage"))

plot_long |> 
  ggplot(aes(x = endpoint, y = value)) +
  geom_segment(aes(xend = endpoint, yend=0), color = "grey", linewidth = 2) +
  geom_point(size=4, color="black") +
  facet_grid(rows = vars(scenario), cols = vars(variable)) +
  coord_flip() +
  theme_bw() +
  xlab("Endpoint ID")

ggsave(filename = glue::glue("plots/lollipop_cf_priors_{target}_data_scenario_1_{model_type}.png"),
       height = 20, width = 20, dpi = 640, units = "cm")

####################
# pair scatter plot
##TODO:
# load(glue::glue("data/summary_data_cluster.RData"))
# summary_dat_hier <- summary_dat
# load(glue::glue("data/summary_data_cluster_separate.RData"))
# summary_dat_sep <- summary_dat



#########
# tables
#########

#' take average across curves
#'
scenario_mean_table <- function(data_list) {
  
  result <- data.frame(Scenario = character(),
                       Coverage = numeric(),
                       Bias = numeric(),
                       RB = numeric(),
                       EmpiricalSE = numeric(),
                       stringsAsFactors = FALSE)
  
  for (i in seq_along(data_list)) {
    temp_data <- data_list[[i]]
    coverage <- mean(temp_data[, "coverage"])
    abs_bias_mean <- mean(abs(temp_data[, "bias"]))
    relative_bias <- mean(temp_data[, "relative_bias"])
    empirical_se <- mean(temp_data[, "empirical_se"])
    
    result <- rbind(result,
                    data.frame(Scenario = i,
                               Coverage = coverage,
                               Bias = abs_bias_mean,
                               RB = relative_bias,
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
  select(data_id, n_endpoints, nsample, prop_censoring, sigma_true) |>
  # filter(data_id %in% 1:16) |>
  kable(format = "latex", booktabs = TRUE)

output_table <-
  kable(combined_tab, format = "latex", booktabs = TRUE) |> 
  add_header_above(c(" ", "RMST" = 4, "Cure Fraction" = 4, "Median" = 4)) 
# kable_styling(latex_options = c("striped", "hold_position"))

# save
write.csv(rmst_tab, file = glue::glue("output_data/rmst_table_{model_type}.csv"), row.names = FALSE)
write.csv(cf_tab, file = glue::glue("output_data/cf_table_{model_type}.csv"), row.names = FALSE)
write.csv(median_tab, file = glue::glue("output_data/median_table_{model_type}.csv"), row.names = FALSE)

##############################
# example survival plots for a single run


