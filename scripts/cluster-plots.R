# cluster plots script


# read in both model types so can compare directly

metrics <- list()
summary_list <- list()

for (model_type in c("separate", "hierarchical")) {
  load(glue::glue("data/performance_measures_cluster_{model_type}.RData"))
  load(glue::glue("data/summary_data_cluster_{model_type}.RData"))
  
  metrics[[model_type]] <- pm
  summary_list[[model_type]] <- summary_dat
}

##########################
# histograms of theta_hat
# mean estimate
##########################

# target <- "rmst"
# target <- "median"
target <- "cf"

tx <- 1
# tx <- 2

# extract target data
hist_dat <- list()
for (i in names(summary_list)) {
  hist_dat[[i]] <- 
    summary_list[[i]] |> 
    map(target) |> 
    map(~ .x[[tx]])  # filter by treatment
}

n_scenarios <- length(hist_dat[[1]])

x_min <- min(map_dbl(hist_dat$separate, ~ min(.x$theta_true)),
             map_dbl(hist_dat$hierarchical, ~ min(.x$theta_true)))
x_max <- max(map_dbl(hist_dat$separate, ~ max(.x$theta_true)),
             map_dbl(hist_dat$hierarchical, ~ max(.x$theta_true)))

plot_list <- list()

for (i in seq_len(n_scenarios)) {
  # combine model data into a single data frame
  theta_hat_sep <- hist_dat$separate[[i]]$theta_hat
  theta_hat_hier <- hist_dat$hierarchical[[i]]$theta_hat
  
  df <-
    data.frame(
      value = theta_hat_hier,
      group = "hier") |> 
    rbind(data.frame(
      value = theta_hat_sep,
      group = "sep")) |> 
    mutate(group = factor(group, levels = c("hier", "sep")))
  
  plot_list[[i]] <-
    df |> 
    ggplot(aes(x = value, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
    # mean true values
    # geom_vline(xintercept = mean(hist_dat$separate[[i]]$theta_true), color = "blue") +
    # geom_vline(xintercept = mean(hist_dat$hierarchical[[i]]$theta_true), color = "red") +
    ggtitle(glue::glue("Scenario {i}")) +
    xlab(target) +
    ylab("Frequency") +
    xlim(x_min, x_max) +
    theme_bw()
}

patchwork::wrap_plots(plot_list, ncol = 4)

