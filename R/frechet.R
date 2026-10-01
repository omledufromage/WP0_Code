# Frechet regression: Y | X ~ Frechet(shape = alpha,
# scale = exp(beta0 + beta1 * X)), fitted by conditional MLE (CMLE).
# Source from any script with: source("R/frechet.R")
# Requires the evd package for data generation.


##############################################################################
# Theoretical Fisher information

# Per-observation Fisher information for theta = (alpha, beta0, beta1).
#
# alpha: shape parameter.
# E, Var: expectation and variance of the covariate X.
fisher_info_frechet <- function(alpha, E, Var) {

  # Euler-Mascheroni constant
  gam <- -digamma(1)

  matrix(
    c(
      1 / alpha^2 * (pi^2 / 6 + (1 - gam)^2),
      (1 - gam),
      (1 - gam)*E,

      1 - gam,
      alpha^2,
      alpha^2*E,

      (1-gam)*E,
      alpha^2*E,
      alpha^2*(Var + E^2)
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("alpha", "beta0", "beta1"),
      c("alpha", "beta0", "beta1")
    )
  )
}


##############################################################################
# Log-likelihood

# Conditional log-likelihood of eta = (alpha, beta0, beta1) given data y, x.
# Always returns a finite value so it can be used with L-BFGS-B.
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
# Data generation and estimation

# Draw n observations from the model.
#
# theta:    c(alpha, beta0, beta1).
# cov_dist: name of the covariate's random number generator, e.g. "rt".
# cov_args: list of arguments for cov_dist, without n.
#
# Returns list(x, y).
sim_frechet_reg <- function(n, theta, cov_dist, cov_args) {

  args <- cov_args
  args$n <- n
  X <- do.call(cov_dist, args)

  sigma <- exp(theta[2] + theta[3] * X)

  Y <- evd::rfrechet(
    n,
    loc = 0,
    scale = sigma,
    shape = theta[1]
  )

  list(x = X, y = Y)
}

# CMLE of (alpha, beta0, beta1) by L-BFGS-B, with alpha >= alpha_lower.
#
# Returns c(alpha_hat, beta0_hat, beta1_hat, convergence), where convergence
# is optim()'s code (0 = success).
fit_frechet_cmle <- function(y, x, alpha_lower = 1e-6) {

  fit <- optim(
    par = c(
      alpha = 1,
      beta0 = mean(log(y)),
      beta1 = 0
    ),
    fn = ll_frechet,
    y = y,
    x = x,
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
}


##############################################################################
# Monte Carlo diagnostics
#
# In the functions below, est is a matrix (or data.frame) of estimates with
# one row per replication and one column per parameter, theta is the vector
# of true values, n is the sample size and fi the per-observation Fisher
# information.

# Bias and RMSE of each column of est.
bias_rmse <- function(est, theta) {

  est <- as.matrix(est)

  errors <- sweep(est, 2, theta, "-")

  rbind(
    bias = colMeans(errors),
    rmse = sqrt(colMeans(errors^2))
  )
}

# sqrt(n) * (est - theta) divided by the asymptotic standard deviations, so
# each column should be approximately N(0, 1).
standardize_estimates <- function(est, theta, n, fi) {

  z <- sweep(
    as.matrix(est),
    2,
    theta,
    "-"
  )

  z <- sqrt(n) * z

  asymptotic_sd <- sqrt(diag(solve(fi)))

  z <- sweep(
    z,
    2,
    asymptotic_sd,
    "/"
  )

  as.data.frame(z)
}

# Squared Mahalanobis distances n * (est - theta)' fi (est - theta), one per
# replication. Approximately chi-squared with ncol(est) degrees of freedom.
mahalanobis_d2 <- function(est, theta, n, fi) {

  errors <- sweep(
    as.matrix(est),
    2,
    theta,
    "-"
  )

  apply(errors, 1, function(e) {
    n * drop(t(e) %*% fi %*% e)
  })
}

# QQ plot of squared Mahalanobis distances D2 against chi-squared quantiles,
# with both axes on [0, lim]. The title reports how many points lie above lim.
mahalanobis_qq <- function(D2, n, lim, df = 3) {

  theoretical <- qchisq(
    ppoints(length(D2)),
    df = df
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
