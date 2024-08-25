# cure fraction prior plots
# for robustness analysis inputs


save_plot <- FALSE

#
sample_zero_truncated_cauchy <- function(n = 1000,
                                         location = 1,
                                         scale = 1) {
  samples <- numeric(n)
  i <- 1
  while (i <= n) {
    sample <- rcauchy(n=1, location, scale)
    if (sample > 0) {
      samples[i] <- sample
      i <- i + 1
    }
  }
  samples
}

# samples <- sample_zero_truncated_cauchy(n = 10000, location = 0.4, scale = 2.4)
samples <- sample_zero_truncated_cauchy(n = 10000, location = 10, scale = 1)
hist(samples[samples < 5000], breaks=20000, xlim = c(0,30), main = "", xlab = "")

summary(samples)
quantile(samples, probs = c(0.025, 0.5, 0.975))

density_est <- density(samples[samples < 5000], n = 10000)
plot(density_est, xlim = c(0,30), main = "", xlab = "$\\sigma$")

#
invlogit <- function(x) exp(x)/(1 + exp(x))

# global cure fraction mean prior
# logistic scale
n <- 100000

params <- 
  list(mean = 
         list(
           list(-1.39, 0.2),
           list(-0.3, 0.2),
           list(-0.3, 0.2),
           list(-3, 0.2),
           list(0, 1)),
       sd = 
         list(
           list(0.4, 2.5),
           list(0, 1/10),
           list(1, 1/3),
           list(0, 1/10),
           list(1, 1)))

titles <- 
  list("Base case",
       "Optimistic informative",
       "Optimistic weakly informative",
       "Pessimistic informative",
       "Weakly informative")

save(params, file = "data/cure_fraction_hyperparameters.RData")

# i <- 1
if (save_plot) {
  png(filename = "plots/cure_fraction_prior_densities.png",
      width = 20, height = 20, units = "cm", res = 640)
}

par(mfrow = c(2,3))

for (i in 1:length(titles)) {
  cf_mean <- do.call(rnorm,
                     c(n, params$mean[[i]]))
  cf_sd <- do.call(sample_zero_truncated_cauchy,
                   c(n, params$sd[[i]]))
  
  cf_endpoint <- rnorm(n, mean = cf_mean, sd = cf_sd)
  cf <- invlogit(cf_endpoint)
  
  # endpoint cure fraction plot
  density_est <- density(cf, na.rm = TRUE, from = 0, to = 1, bw = 0.1)
  
  plot_title <- bquote(
    atop(
    .(gt::vec_fmt_roman(i, case = "lower")) * ") " * .(titles[[i]]), 
      "logit(" * pi * ") ~ N(-1.39, " * 0.2^2 * "), " * sigma * " ~ tCauchy(0.4, " * 2.5^2 * ")"))

  plot(density_est, col = "blue", lwd = 2,
       main = plot_title,
       xlab = "Probabiltity",
       cex.lab=1.5,  cex.main=1.5)
  abline(v = 0.2, col = "red", lwd = 2)
}

if (save_plot) dev.off()
