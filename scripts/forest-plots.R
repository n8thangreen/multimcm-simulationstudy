## based on ChatGPT code

library(rstan)
library(brms)
library(broom)
library(dplyr)
library(ggplot2)


load(file = "data/stan_out.RData")

# Extract parameter summaries
summary_fit1 <- 
  as.data.frame(summary(stan_out[[1]]$output)$summary) |> 
  tibble::rownames_to_column() |>
  filter(grepl("^cf", rowname)) |>
  filter(grepl("\\[1", rowname)) |> 
  mutate(model = 1)

summary_fit2 <- 
  as.data.frame(summary(stan_out[[2]]$output)$summary) |> 
  tibble::rownames_to_column() |> 
  filter(grepl("^cf", rowname)) |> 
  filter(grepl("\\[1", rowname)) |> 
  mutate(model = 2)

# summary_fit3 <- as.data.frame(summary(stan_out[[2]]$output)$summary) |> tibble::rownames_to_column() |> filter(grepl("^cf", rowname)) |> filter(grepl("\\[1", rowname)) |> mutate(model = 3)

# Prepare data for plotting
parameters <- 
  do.call(rbind, list(summary_fit1, summary_fit2)) |> 
  mutate(model = as.factor(model))

# Create forest plot
ggplot(parameters, aes(x = rowname, y = mean, color = model)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = 0.5)) +
  theme_minimal() +
  coord_flip() +
  labs(# title = "Forest Plot of Model Estimates",
    x = "Parameter", y = "Probability") +
  ylim(0, 0.5)

##TODO: add true value to plot
