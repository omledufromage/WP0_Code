##############################################################################
# CMLE for Frechet regression: Y | X ~ Frechet(shape = alpha,
# scale = exp(beta0 + beta1 * X)).
# Run from the project root (opening Code.Rproj does this for you):
#   Rscript simulations/01_frechet_cmle.R

source("R/utils.R")

##############################################################################
# Parameters

# TRUE: run the simulation and save it. FALSE: load the latest saved run.
run_sim <- TRUE

sample_sizes <- c(50, 250, 500, 1000, 2000)
# sample_sizes <- c(50, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000)
M <- 5000

alpha_true <- exp(1)
beta0_true <- pi
beta1_true <- sqrt(2)

theta_true <- c(alpha_true, beta0_true, beta1_true)

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
E        <- covariates[[covariate]]$E
Var      <- covariates[[covariate]]$Var

study <- paste0("frechet_cmle_", covariate)

# Sample sizes shown in the QQ plots (must be in sample_sizes)
qq_n  <- 500
mah_n <- c(50, 500, 1000, 2000)
stopifnot(all(c(qq_n, mah_n) %in% sample_sizes))

# Euler-Mascheroni constant
gam <- -digamma(1)

# Theoretical Fisher Information
fi <- matrix(
  c(
    1 / alpha_true^2 * (pi^2 / 6 + (1 - gam)^2),
    (1 - gam),
    (1 - gam)*E,

    1 - gam,
    alpha_true^2,
    alpha_true^2*E,

    (1-gam)*E,
    alpha_true^2*E,
    alpha_true^2*(Var + E^2)
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c("alpha", "beta0", "beta1"),
    c("alpha", "beta0", "beta1")
  )
)


##############################################################################
# Log-likelihood

ll_frechet <- function(eta, y, x) {

  alpha <- eta[1]
  beta0 <- eta[2]
  beta1 <- eta[3]

  if (alpha <= 0){
    s <- (-Inf)
  } else {
    d <- beta0 + beta1 * x - log(y)

    s <- sum(
      log(alpha) - log(y) +
        alpha * d -
        exp(alpha * d)
    )
  }

  # L-BFGS-B needs finite values. NaN is treated like -Inf.
  if (!is.finite(s)) {
    if (is.nan(s) || s < 0) {
      cat("negative infinity happened\n")
      s <- -10^308
    } else {
      cat("INFINITY happened\n")
      s <- 10^308
    }
  }

  s
}


##############################################################################
# Monte Carlo simulation

simulate <- function() {

  set.seed(02072024)

  results <- lapply(sample_sizes, function(n) {

    cat("Running n =", n, "\n")

    experiment <- replicate(M, {

      args <- cov_args
      args$n <- n
      X <- do.call(cov_dist, args)

      sigma <- exp(beta0_true + beta1_true * X)

      Y <- evd::rfrechet(
        n,
        loc = 0,
        scale = sigma,
        shape = alpha_true
      )

      # CMLE
      fit <- optim(
        par = c(
          alpha = 1,
          beta0 = mean(log(Y)),
          beta1 = 0
        ),
        fn = ll_frechet,
        y = Y,
        x = X,
        method = "L-BFGS-B",
        lower = c(alpha_lower, -Inf, -Inf),
        upper = c(Inf, Inf, Inf),
        control = list(fnscale = -1)
      )

      c(
        alpha_hat = unname(fit$par[1]),
        beta0_hat = unname(fit$par[2]),
        beta1_hat = unname(fit$par[3]),
        convergence = fit$convergence
      )

    })

    # replicate() gives parameters in rows, simulations in columns
    experiment <- t(experiment)
    experiment <- as.data.frame(experiment)

    # Keep track of optimization failures, including fits stuck on the
    # alpha lower bound (optim reports those as converged)
    at_bound <- experiment$alpha_hat <= 2 * alpha_lower
    ok <- experiment$convergence == 0 & !at_bound

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
      n_at_bound = sum(at_bound)
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
      theta_true = theta_true,
      covariate = covariates[[covariate]],
      alpha_lower = alpha_lower,
      seed = 02072024
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
  n_at_bound = sapply(results, function(res) res$n_at_bound)
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

  estimates <- as.matrix(res$estimates)

  bias <- colMeans(estimates) - true_values

  rmse <- sqrt(
    colMeans(
      sweep(estimates, 2, true_values, "-")^2
    )
  )

  c(
    n = res$n,

    bias_alpha = unname(bias["alpha_hat"]),
    bias_beta0 = unname(bias["beta0_hat"]),
    bias_beta1 = unname(bias["beta1_hat"]),

    rmse_alpha = unname(rmse["alpha_hat"]),
    rmse_beta0 = unname(rmse["beta0_hat"]),
    rmse_beta1 = unname(rmse["beta1_hat"])
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

fi_inv <- solve(fi)

standardized <- lapply(results, function(res) {

  est <- as.matrix(res$estimates)
  n <- res$n

  z <- sweep(
    est,
    2,
    true_values,
    "-"
  )

  z <- sqrt(n) * z

  asymptotic_sd <- sqrt(diag(fi_inv))

  z <- sweep(
    z,
    2,
    asymptotic_sd,
    "/"
  )

  as.data.frame(z)
})

names(standardized) <- names(results)

z_qq <- standardized[[paste0("n_", qq_n)]]

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

mahalanobis_d2 <- function(res, fi, theta_true) {

  est <- as.matrix(res$estimates)

  errors <- sweep(
    est,
    2,
    theta_true,
    "-"
  )

  apply(errors, 1, function(e) {
    res$n * drop(t(e) %*% fi %*% e)
  })
}

mahalanobis_qq <- function(D2, n, lim) {

  theoretical <- qchisq(
    ppoints(length(D2)),
    df = 3
  )

  # Points beyond lim are not drawn, so report how many there are
  n_out <- sum(D2 > lim)

  plot(
    theoretical,
    sort(D2),
    xlab = "",
    ylab = "",
    main = paste0("n = ", n, if (n_out > 0) paste0(" (", n_out, " above ", round(lim), ")")),
    pch = 16,
    cex = 0.5,
    col = adjustcolor("blue", alpha.f = 0.3),
    xlim = c(0, lim),
    ylim = c(0, lim)
  )

  abline(0, 1, lty = 2)
}

D2_list <- lapply(
  results[paste0("n_", mah_n)],
  mahalanobis_d2,
  fi = fi,
  theta_true = theta_true
)

# Common axis limit for all panels: the 99.5% quantile of D2, but at least
# the largest theoretical quantile
lim <- max(
  qchisq(ppoints(M)[M], df = 3),
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
