
#
sample_zero_truncated_cauchy <- function(n = 1000, location = 1, scale = 1) {
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

samples <- sample_zero_truncated_cauchy(n = 1000000, location = 0.4, scale = 2.4)
hist(samples, breaks=100000, xlim = c(0,50))

summary(samples)
quantile(samples, probs = c(0.025, 0.5, 0.975))

density_est <- density(samples)
plot(density_est)

#
invlogit <- function(x) exp(x)/(1 + exp(x))

# global cure fraction mean prior
# logistic scale
n <- 10000

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

png(filename = "plots/cure_fraction_prior_densities.png",
    width = 20, height = 20, units = "cm", res = 640)

par(mfrow = c(3,2))

for (i in 1:5) {
  cf_mean <- do.call(rnorm, c(n, params$mean[[i]]))
  cf_sd <- do.call(sample_zero_truncated_cauchy, c(n, params$sd[[i]]))
  
  cf_global <- rnorm(n, cf_mean, cf_sd)
  cf <- invlogit(cf_global)
  
  # global cure fraction plot
  # hist(cf, breaks = 100, freq = FALSE, main = titles[[i]])
  density_est <- density(cf, na.rm = TRUE, from = 0, to = 1, bw = 0.1)
  plot(density_est, col = "blue", lwd = 2,
       main = paste0(gt::vec_fmt_roman(i, case = "lower"), ") ", titles[[i]]),
       xlab = "Probability")
  abline(v = 0.2, col = "red", lwd = 2)
}

dev.off()
