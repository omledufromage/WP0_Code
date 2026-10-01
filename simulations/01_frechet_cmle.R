##############################################################################
# CMLE for Frechet regression: Y | X ~ Frechet(shape = alpha,
# scale = exp(beta0 + beta1 * X)).
# Model functions are in R/frechet.R.
# Run from the project root (opening Code.Rproj does this for you):
#   Rscript simulations/01_frechet_cmle.R

source("R/utils.R")
source("R/frechet.R")

##############################################################################
# Parameters

# TRUE: run the simulation and save it. FALSE: load the latest saved run.
run_sim <- TRUE

sample_sizes <- c(50, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000)
# sample_sizes <- c(50, 250, 500, 1000, 2000)
M <- 5000
seed <- 02072024

alpha_true <- exp(1)
beta0_true <- pi
beta1_true <- sqrt(2)

true_values <- c(
  alpha_hat = alpha_true,
  beta0_hat = beta0_true,
  beta1_hat = beta1_true
)

# Lower bound on alpha in optim(). Fits that end on it count as failures.
alpha_lower <- 1e-6

# Covariate distributions, with the expectation and variance of the
# random variable with distribution dist and arguments args.
covariates <- list(
  exp4 = list(dist = "rexp", args = list(rate = 4), E = 1/4, Var = 1/16),
  t6   = list(dist = "rt",   args = list(df = 6),   E = 0,   Var = 6/(6 - 2))
)

covariate <- "t6"

cov_dist <- covariates[[covariate]]$dist
cov_args <- covariates[[covariate]]$args

study <- paste0("frechet_cmle_", covariate)

# Sample sizes shown in the QQ plots (must be in sample_sizes)
qq_n  <- 500
mah_n <- c(50, 500, 1000, 2000)
stopifnot(all(c(qq_n, mah_n) %in% sample_sizes))

# Theoretical Fisher Information
fi <- fisher_info_frechet(
  alpha_true,
  E = covariates[[covariate]]$E,
  Var = covariates[[covariate]]$Var
)


##############################################################################
# Monte Carlo simulation

simulate <- function() {

  set.seed(seed)

  results <- lapply(sample_sizes, function(n) {

    cat("Running n =", n, "\n")

    experiment <- replicate(M, {
      dat <- sim_frechet_reg(n, true_values, cov_dist, cov_args)
      fit_frechet_cmle(dat$y, dat$x, alpha_lower)
    })

    # replicate() gives parameters in rows, simulations in columns
    experiment <- t(experiment)
    experiment <- as.data.frame(experiment)

    # Keep track of optimization failures, including fits stuck on the
    # alpha lower bound (optim reports those as converged). Fits where
    # optim stopped with an error have convergence = -1.
    at_bound <- !is.na(experiment$alpha_hat) &
      experiment$alpha_hat <= 2 * alpha_lower
    ok <- experiment$convergence == 0 & !at_bound
    n_errors <- sum(experiment$convergence == -1)

    convergence_rate <- mean(ok)

    # Parameter estimates only
    estimates <- experiment[ok, c("alpha_hat", "beta0_hat", "beta1_hat")]

    mu <- colMeans(estimates)
    s <- apply(estimates, 2, sd)

    # Monte Carlo estimate of per-observation Fisher information
    fi_mc <- solve(cov(estimates)) / n

    list(
      n = n,
      estimates = estimates,
      mu = mu,
      sd = s,
      fi_mc = fi_mc,
      convergence_rate = convergence_rate,
      n_at_bound = sum(at_bound),
      n_errors = n_errors
    )
  })

  names(results) <- paste0("n_", sample_sizes)
  results
}

if (run_sim) {
  results <- simulate()
  save_results(
    results,
    study,
    settings = list(
      M = M,
      sample_sizes = sample_sizes,
      theta_true = true_values,
      covariate = covariates[[covariate]],
      alpha_lower = alpha_lower,
      seed = seed
    )
  )
} else {
  saved <- load_latest_results(study)
  results <- saved$results
  sample_sizes <- saved$settings$sample_sizes
  stopifnot(all(c(qq_n, mah_n) %in% sample_sizes))
}

##############################################################################
# Summary of parameter estimates

mean_table <- do.call(
  rbind,
  lapply(results, function(res) res$mu)
)

sd_table <- do.call(
  rbind,
  lapply(results, function(res) res$sd)
)

rownames(mean_table) <- sample_sizes
rownames(sd_table) <- sample_sizes

mean_table
sd_table

bias_table <- sweep(
  mean_table,
  2,
  true_values,
  "-"
)

bias_table

################################################################################
# Checking convergence

convergence_table <- rbind(
  convergence_rate = sapply(results, function(res) res$convergence_rate),
  n_at_bound = sapply(results, function(res) res$n_at_bound),
  n_errors = sapply(results, function(res) res$n_errors)
)

convergence_table

################################################################################
# Compare MC simulation to theoretical Fisher Information

for (i in seq_along(results)) {

  cat("\n----------------------------------\n")
  cat("n =", results[[i]]$n, "\n")
  cat("----------------------------------\n")

  print(results[[i]]$fi_mc)
}

fi_error <- lapply(results, function(res) {
  res$fi_mc - fi
})

fi_error

################################################################################
# Bias and RMSE

summary_table <- do.call(rbind, lapply(results, function(res) {

  br <- bias_rmse(res$estimates, true_values)

  c(
    n = res$n,

    bias_alpha = unname(br["bias", "alpha_hat"]),
    bias_beta0 = unname(br["bias", "beta0_hat"]),
    bias_beta1 = unname(br["bias", "beta1_hat"]),

    rmse_alpha = unname(br["rmse", "alpha_hat"]),
    rmse_beta0 = unname(br["rmse", "beta0_hat"]),
    rmse_beta1 = unname(br["rmse", "beta1_hat"])
  )
}))

summary_table

cols <- 1:3

param_labels <- c(
  expression(hat(alpha)),
  expression(hat(beta)[0]),
  expression(hat(beta)[1])
)

# RMSE
matplot(
  summary_table[, "n"],
  summary_table[, c("rmse_alpha", "rmse_beta0", "rmse_beta1")],
  type = "b",
  pch = 1,
  lty = 1,
  col = cols,
  #log = "xy",
  xlab = "Sample size n",
  ylab = "RMSE",
  main = "RMSE of the CMLE"
)

legend(
  "topright",
  legend = param_labels,
  col = cols,
  lty = 1,
  pch = 1,
  bty = "n"
)

# Bias
matplot(
  summary_table[, "n"],
  summary_table[, c("bias_alpha", "bias_beta0", "bias_beta1")],
  type = "b",
  pch = 1,
  lty = 2,
  col = cols,
  xlab = "Sample size n",
  ylab = "Bias",
  main = "Bias of the CMLE"
)

abline(h = 0, lty = 3)

legend(
  "topright",
  legend = param_labels,
  col = cols,
  lty = 2,
  pch = 1,
  bty = "n"
)

################################################################################
# Marginal QQ plots of sqrt(n) * (theta_hat - theta) / asymptotic sd

res_qq <- results[[paste0("n_", qq_n)]]
z_qq <- standardize_estimates(res_qq$estimates, true_values, qq_n, fi)

old_par <- par(mfrow = c(1, 3))

for (j in 1:3) {

  qqnorm(
    z_qq[[j]],
    main = paste0(names(z_qq)[j], ", n = ", qq_n),
    xlab = "Theoretical normal quantiles",
    ylab = "Empirical quantiles"
  )

  qqline(
    z_qq[[j]]
  )
}

par(old_par)

################################################################################
# Mahalanobis distances

D2_list <- lapply(results[paste0("n_", mah_n)], function(res) {
  mahalanobis_d2(res$estimates, true_values, res$n, fi)
})

# Common axis limit for all panels: the 99.5% quantile of D2, but at least
# the largest theoretical quantile
max_reps <- max(lengths(D2_list))
lim <- max(
  qchisq(ppoints(max_reps)[max_reps], df = 3),
  quantile(unlist(D2_list), 0.995)
)

old_par <- par(
  mfrow = c(2, 2),
  mar = c(2, 2, 2, 1),
  oma = c(4, 4, 1, 1)
)

for (i in seq_along(mah_n)) {
  mahalanobis_qq(D2_list[[i]], mah_n[i], lim)
}

mtext(
  expression("Theoretical " * chi[3]^2 * " quantiles"),
  side = 1,
  outer = TRUE,
  line = 2
)

mtext(
  expression("Empirical squared Mahalanobis distance " * D^2),
  side = 2,
  outer = TRUE,
  line = 2
)

par(old_par)
