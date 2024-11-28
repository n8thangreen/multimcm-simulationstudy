# tables of performance measures output


#' average across curves
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


library(kableExtra)
library(purrr)
library(dplyr)

# load data

model_type <- "separate"

load(glue::glue("data/performance_measures_cluster_{model_type}.RData"))
load(glue::glue("data/summary_data_cluster_{model_type}.RData"))

pm_sep <- pm
summary_dat_sep <- summary_dat

model_type <- "hierarchical"

load(glue::glue("data/performance_measures_cluster_{model_type}.RData"))
load(glue::glue("data/summary_data_cluster_{model_type}.RData"))

pm_hier <- pm
summary_dat_hier <- summary_dat

## both model types in a single table

rmst_tab_sep <- 
  map(pm_sep, "rmst") |> 
  scenario_mean_table()

cf_tab_sep <- 
  map(pm_sep, "cf") |> 
  scenario_mean_table()

rmst_tab_hier <- 
  map(pm_hier, "rmst") |> 
  scenario_mean_table()

cf_tab_hier <- 
  map(pm_hier, "cf") |> 
  scenario_mean_table()

combined_cf_tab <- 
  merge(cf_tab_sep, cf_tab_hier,
        by = "Scenario",
        suffixes = c(".sep", ".hier")) |> 
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(order(colnames(.))) |> 
  select(Scenario, everything())

write.csv(combined_cf_tab, file = glue::glue("output_data/cf_table.csv"), row.names = FALSE)

combined_rmst_tab <- 
  merge(rmst_tab_sep, rmst_tab_hier,
        by = "Scenario",
        suffixes = c(".sep", ".hier")) |> 
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(order(colnames(.))) |> 
  select(Scenario, everything())

write.csv(combined_rmst_tab, file = glue::glue("output_data/rmst_table.csv"), row.names = FALSE)

## single table

combined_tab_all <- 
  merge(combined_rmst_tab, combined_cf_tab, by = "Scenario")

write.csv(combined_tab_all, file = glue::glue("output_data/combined_performance_table.csv"), row.names = FALSE)

## single table per model type

rmst_tab <- 
  map(pm, "rmst") |> 
  scenario_mean_table()

cf_tab <- 
  map(pm, "cf") |> 
  scenario_mean_table()

combined_tab <- 
  merge(rmst_tab, cf_tab, by = "Scenario") |> 
  mutate(across(where(is.numeric), ~ round(., 2)))

#############
# formatting

library(stringr)

# xtable::xtable(combined_tab, digits = 3)

##TODO:
# input_table <-
#   scenario_data |> 
#   select(data_id, n_endpoints, nsample, prop_censoring, sigma_true) |>
#   # filter(data_id %in% 1:16) |>
#   kable(format = "latex", booktabs = TRUE)

output_table <-
  kable(combined_tab, format = "latex", booktabs = TRUE) |> 
  add_header_above(c(" ", "RMST" = 4, "Cure Fraction" = 4)) 
# kable_styling(latex_options = c("striped", "hold_position"))

output_table_all <-
  combined_tab_all |> 
  kable(format = "latex", booktabs = TRUE,
        col.names = c("Scenario", rep(c("hier", "sep"), 8))) |> 
  add_header_above(c(" ",
                     "Bias" = 2, "Coverage" = 2, "empSE" = 2, "RB" = 2,
                     "Bias" = 2, "Coverage" = 2, "emp SE" = 2, "RB" = 2)) |>
  add_header_above(c(" ", "RMST" = 8, "Cure Fraction" = 8))
  

# save
write.csv(rmst_tab, file = glue::glue("output_data/rmst_table_{model_type}.csv"), row.names = FALSE)
write.csv(cf_tab, file = glue::glue("output_data/cf_table_{model_type}.csv"), row.names = FALSE)
write.csv(median_tab, file = glue::glue("output_data/median_table_{model_type}.csv"), row.names = FALSE)

