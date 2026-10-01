# Frechet regression: Y | X ~ Frechet(shape = alpha,
# scale = exp(beta0 + beta1 * X1 + ... + betap * Xp)), fitted by conditional
# MLE (CMLE). The parameter vector is theta = (alpha, beta0, ..., betap).
# Source from any script with: source("R/frechet.R")
# Requires the evd package for data generation.


# Parameter names for a model with p covariates.
frechet_par_names <- function(p) {
  c("alpha", paste0("beta", 0:p))
}


##############################################################################
# Theoretical Fisher information

# Per-observation Fisher information for theta = (alpha, beta0, ..., betap).
#
# alpha: shape parameter.
# E:     vector of the covariates' expectations (length p).
# Var:   vector of the covariates' variances (independent covariates), or
#        their p x p covariance matrix.
fisher_info_frechet <- function(alpha, E, Var) {

  # Euler-Mascheroni constant
  gam <- -digamma(1)

  p <- length(E)
  Sigma <- if (is.matrix(Var)) Var else diag(Var, nrow = p)

  # E[z] and E[z z'] for z = (1, X)
  Ez <- c(1, E)
  Ezz <- rbind(
    c(1, E),
    cbind(E, Sigma + E %*% t(E))
  )

  fi <- matrix(0, p + 2, p + 2)
  fi[1, 1] <- 1 / alpha^2 * (pi^2 / 6 + (1 - gam)^2)
  fi[1, -1] <- (1 - gam) * Ez
  fi[-1, 1] <- (1 - gam) * Ez
  fi[-1, -1] <- alpha^2 * Ezz

  dimnames(fi) <- list(frechet_par_names(p), frechet_par_names(p))
  fi
}


##############################################################################
# Log-likelihood

# Conditional log-likelihood of eta = (alpha, beta0, ..., betap) given data
# y and the n x p covariate matrix x (a vector is treated as p = 1).
# Always returns a finite value so it can be used with L-BFGS-B.
ll_frechet <- function(eta, y, x) {

  x <- as.matrix(x)

  alpha <- eta[1]
  beta0 <- eta[2]
  beta <- eta[-(1:2)]

  if (alpha <= 0){
    s <- (-Inf)
  } else {
    d <- drop(beta0 + x %*% beta) - log(y)

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

# Draw n observations from the model with independent covariates.
#
# theta: c(alpha, beta0, beta1, ..., betap).
# covs:  list of p covariate specifications, each a list with
#        dist (name of a random number generator, e.g. "rt") and
#        args (list of its arguments, without n).
#
# Returns list(x, y), where x is the n x p covariate matrix.
sim_frechet_reg <- function(n, theta, covs) {

  X <- sapply(covs, function(cv) {
    args <- cv$args
    args$n <- n
    do.call(cv$dist, args)
  })
  X <- matrix(X, nrow = n, dimnames = list(NULL, paste0("x", seq_along(covs))))

  sigma <- exp(drop(theta[2] + X %*% theta[-(1:2)]))

  Y <- evd::rfrechet(
    n,
    loc = 0,
    scale = sigma,
    shape = theta[1]
  )

  list(x = X, y = Y)
}

# Moment-based starting values from a least-squares fit of log(y) on x.
# Under the model, log(Y) = beta0 + X'beta + G / alpha with G standard
# Gumbel (mean gam, variance pi^2 / 6).
start_frechet <- function(y, x, alpha_lower = 1e-6) {

  gam <- -digamma(1)

  ols <- lm.fit(cbind(1, x), log(y))
  alpha0 <- max(pi / (sqrt(6) * sd(ols$residuals)), 2 * alpha_lower)
  beta_ols <- unname(ols$coefficients)

  c(alpha0, beta_ols[1] - gam / alpha0, beta_ols[-1])
}

# CMLE of (alpha, beta0, ..., betap) by L-BFGS-B, with alpha >= alpha_lower,
# started from start_frechet().
#
# Returns c(alpha_hat, beta0_hat, ..., betap_hat, convergence), where
# convergence is optim()'s code (0 = success), or -1 if optim() stopped with
# an error (the estimates are then NA).
fit_frechet_cmle <- function(y, x, alpha_lower = 1e-6) {

  x <- as.matrix(x)
  p <- ncol(x)

  fit <- tryCatch(
    optim(
      par = start_frechet(y, x, alpha_lower),
      fn = ll_frechet,
      y = y,
      x = x,
      method = "L-BFGS-B",
      lower = c(alpha_lower, rep(-Inf, p + 1)),
      upper = rep(Inf, p + 2),
      control = list(fnscale = -1)
    ),
    error = function(e) list(par = rep(NA_real_, p + 2), convergence = -1)
  )

  c(
    setNames(unname(fit$par), paste0(frechet_par_names(p), "_hat")),
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

# QQ plot of squared Mahalanobis distances D2 against chi-squared quantiles
# with df degrees of freedom (the number of parameters), with both axes on
# [0, lim]. The title reports how many points lie above lim.
mahalanobis_qq <- function(D2, n, lim, df) {

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
