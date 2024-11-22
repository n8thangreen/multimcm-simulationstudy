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

target <- "rmst"
# target <- "median"
# target <- "cf"

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
             map_dbl(hist_dat$separate, ~ min(.x$theta_hat)),
             map_dbl(hist_dat$hierarchical, ~ min(.x$theta_true)),
             map_dbl(hist_dat$hierarchical, ~ min(.x$theta_hat)))

x_max <- max(map_dbl(hist_dat$separate, ~ max(.x$theta_true)),
             map_dbl(hist_dat$separate, ~ max(.x$theta_hat)),
             map_dbl(hist_dat$hierarchical, ~ max(.x$theta_true)),
             map_dbl(hist_dat$hierarchical, ~ max(.x$theta_hat)))

plot_list <- list()

for (i in seq_len(n_scenarios)) {
  # combine model data into a single data frame
  theta_hat_sep <- hist_dat$separate[[i]]$theta_hat
  theta_hat_hier <- hist_dat$hierarchical[[i]]$theta_hat
  theta_true_sep <- hist_dat$separate[[i]]$theta_true
  theta_true_hier <- hist_dat$hierarchical[[i]]$theta_true
  
  df <-
    data.frame(
      value = theta_hat_hier,
      true = theta_true_hier,
      group = "hier") |> 
    rbind(data.frame(
      value = theta_hat_sep,
      true = theta_true_sep,
      group = "sep")) |> 
    mutate(group = factor(group, levels = c("hier", "sep")))
  
  plot_list[[i]] <-
    df |> 
    ggplot(aes(x = value, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
    geom_density(aes(x = true), col = "darkgreen", linewidth = 1, inherit.aes = FALSE) +
    # mean true values
    # geom_vline(xintercept = mean(hist_dat$separate[[i]]$theta_true), color = "blue") +
    # geom_vline(xintercept = mean(hist_dat$hierarchical[[i]]$theta_true), color = "red") +
    ggtitle(glue::glue("Scenario {i}")) +
    xlab(target) +
    ylab("Frequency") +
    xlim(x_min, x_max) +
    theme_bw() #+
    # theme(legend.position = "none")
}

# Keep legend for the first plot, remove from others
plot_list <- lapply(seq_along(plot_list), function(i) {
  if (i == 1) {
    plot_list[[i]] # Keep the legend for the first plot
  } else {
    plot_list[[i]] + theme(legend.position = "none") # Remove legend
  }
})

patchwork::wrap_plots(plot_list, nrow = 4) +
  patchwork::plot_layout(guides = "collect")

ggsave(filename = glue::glue("plots/theta_hat_hist_{target}.png"),
       height = 20, width = 35, dpi = 640, units = "cm")


#################
# lollipop plots

target <- "rmst"
# target <- "cf"
# target <- "median"

# quo_measure <- quo(bias)
# quo_measure <- quo(relative_bias)
quo_measure <- quo(coverage)
# quo_measure <- quo(empirical_se)

label_data <-
  scenario_data |> 
  mutate(label = glue::glue("{data_id}: e={n_endpoints} n={nsample}")) |> 
  rename(scenario = data_id) |> 
  mutate(label = factor(label, levels = unique(label))) |> 
  select(label, scenario)

# extract target data and combine scenarios
plot_dat_sep <-
  metrics$separate |> 
  map(target) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario") |> 
  mutate(endpoint = factor(endpoint),
         group = "sep") |> 
  left_join(label_data, by = "scenario")

plot_dat_hier <-
  metrics$hierarchical |> 
  map(target) |> 
  map(~as.data.frame(.x)) |> 
  map(~mutate(.x, endpoint = 1:n())) |> 
  list_rbind(names_to = "scenario") |> 
  mutate(endpoint = factor(endpoint),
         group = "hier") |> 
  left_join(label_data, by = "scenario")

plot_dat <- rbind(plot_dat_sep, plot_dat_hier)

# clean measure string
measure <- quo_name(quo_measure)
y_label <- stringr::str_to_sentence(gsub(pattern = "\\_", " ", measure))

plot_dat |> 
  ggplot(aes(x = endpoint, y = !!quo_measure, col = group, fill = group)) +
  # geom_segment(aes(xend = endpoint, yend=0), linewidth = 2, position = position_dodge(width = 0.8)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(size=4, position = position_dodge(width = 0.8)) +
  facet_wrap(vars(label), nrow = 4) +
  # facet_wrap(vars(scenario)) +  # without titles
  coord_flip() +
  theme_minimal() +
  xlab("Endpoint ID") +
  ylab(y_label) #+
  # theme(legend.position = "none")

ggsave(filename = glue::glue("plots/lollipop_{target}_{quo_name(quo_measure)}.png"),
       height = 20, width = 30, dpi = 640, units = "cm", bg = "white")
