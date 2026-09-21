# Data generation for hybrid-control Cox simulations =====================

#' Sample from a multivariate normal via Cholesky factorization
#'
#' @param n Number of samples.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @return An \code{n} x \code{length(mu)} matrix of samples.
#' @keywords internal
#' @noRd
rmvnorm_chol <- function(n, mu, Sigma) {
  p <- length(mu)
  if (p == 0L) return(matrix(numeric(0), nrow = n, ncol = 0))
  L <- tryCatch(
    chol(Sigma),
    error = function(e) {
      warning("Sigma not positive definite; adding small jitter to diagonal")
      chol(Sigma + diag(NUMERIC_ZERO * 100, p))
    }
  )
  Z <- matrix(stats::rnorm(n * p), n, p)
  sweep(Z %*% L, 2, mu, "+")
}

#' Simulate Weibull-proportional-hazards survival times
#'
#' @param n Number of subjects.
#' @param lp Vector of linear predictors of length \code{n}.
#' @param shape Weibull shape parameter (positive).
#' @param lambda Baseline scale parameter (positive).
#' @return Numeric vector of event times.
#' @keywords internal
#' @noRd
sim_weibull_ph <- function(n, lp, shape = 1.2, lambda = 0.02) {
  u <- stats::runif(n)
  (-log(u) / (lambda * exp(lp)))^(1 / shape)
}

#' Bisect for an exponential censoring rate matching a target censoring rate
#'
#' @param T_event Vector of event times.
#' @param target_cens Target censoring proportion.
#' @param max_rate Upper bound for the bisection.
#' @param tol Convergence tolerance on censoring proportion.
#' @param it Maximum bisection iterations.
#' @return The estimated exponential rate.
#' @keywords internal
#' @noRd
pick_censor_rate <- function(T_event, target_cens = 0.2,
                             max_rate = 2, tol = 0.01, it = 30) {
  lo <- 1e-8
  hi <- max_rate
  for (k in 1:it) {
    mid <- (lo + hi) / 2
    C <- stats::rexp(length(T_event), rate = mid)
    cens <- mean(C < T_event)
    if (cens > target_cens) hi <- mid else lo <- mid
    if (abs(cens - target_cens) < tol) break
  }
  (lo + hi) / 2
}

#' Simulate a hybrid-control Cox proportional hazards dataset
#'
#' Generates a randomized trial augmented with an external control arm,
#' under a Cox proportional hazards model
#' \deqn{h(t \mid T,Z,X) = h_0(t)\,\exp(\theta T + \delta Z + \beta^\top X),}
#' where \eqn{T} is the treatment indicator, \eqn{Z} is the external-control
#' indicator, \eqn{X} is a vector of baseline covariates, \eqn{\theta} is
#' the treatment effect, and \eqn{\delta} is the population drift between
#' concurrent and external controls.
#'
#' @param nI1 Number of internal randomized treated subjects.
#' @param nI0 Number of internal randomized concurrent control subjects.
#' @param nE  Number of external control subjects.
#' @param theta0 True treatment effect (log hazard ratio).
#' @param delta0 True population drift (log hazard ratio for external vs.
#'   internal controls).
#' @param p Number of covariates; use 0 for the no-covariate setting.
#' @param beta Numeric vector of covariate coefficients (length \code{p}).
#'   Defaults to \code{rep(0.2, p)}.
#' @param rho Equicorrelation of covariates.
#' @param cov_shift Covariate mean shift in the external control arm
#'   (numeric vector of length \code{p}).
#' @param shape Weibull shape parameter for event-time generation.
#' @param lambda Baseline scale parameter for event-time generation.
#' @param target_cens Target right-censoring proportion.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{data}}{A data frame with \code{time}, \code{status},
#'       \code{T}, \code{Z}, and covariates \code{X1}, ..., \code{Xp}.}
#'     \item{\code{truth}}{The true parameters used to generate the data.}
#'     \item{\code{settings}}{The simulation settings, including the
#'       calibrated censoring rate.}
#'   }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_cox(nI1 = 100, nI0 = 100, nE = 200,
#'                            theta0 = log(0.8), delta0 = 0)
#' head(sim$data)
#'
#' @export
simulate_hybrid_cox <- function(nI1 = 150, nI0 = 150, nE = 300,
                                theta0 = 0, delta0 = 0,
                                p = 5, beta = NULL,
                                rho = 0.0,
                                cov_shift = rep(0, p),
                                shape = 1.2, lambda = 0.02,
                                target_cens = 0.2) {
  stopifnot(
    nI1 > 0, nI0 > 0, nE > 0,
    p >= 0, p == as.integer(p), rho >= -1, rho <= 1,
    length(cov_shift) == p,
    shape > 0, lambda > 0,
    target_cens >= 0, target_cens < 1
  )
  if (nI1 < 10 || nI0 < 10 || nE < 10) {
    warning("Small sample sizes (n < 10) may cause numerical instability")
  }

  if (is.null(beta)) beta <- rep(0.2, p)
  if (length(beta) != p) {
    stop("Length of beta (", length(beta), ") must equal p (", p, ")")
  }

  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  X_I1 <- rmvnorm_chol(nI1, mu = rep(0, p), Sigma = Sigma)
  X_I0 <- rmvnorm_chol(nI0, mu = rep(0, p), Sigma = Sigma)
  X_E  <- rmvnorm_chol(nE,  mu = cov_shift,  Sigma = Sigma)

  X <- rbind(X_I1, X_I0, X_E)
  T <- c(rep(1, nI1), rep(0, nI0), rep(0, nE))
  Z <- c(rep(0, nI1 + nI0), rep(1, nE))

  lp <- as.numeric(theta0 * T + delta0 * Z + X %*% beta)
  T_event <- sim_weibull_ph(n = nI1 + nI0 + nE, lp = lp,
                            shape = shape, lambda = lambda)

  rate_c <- pick_censor_rate(T_event[Z == 0], target_cens = target_cens)
  C <- stats::rexp(length(T_event), rate = rate_c)

  time   <- pmin(T_event, C)
  status <- as.integer(T_event <= C)

  dat <- data.frame(time = time, status = status, T = T, Z = Z)
  for (j in seq_len(p)) dat[[paste0("X", j)]] <- X[, j]

  list(
    data = dat,
    truth = list(theta0 = theta0, delta0 = delta0, beta = beta),
    settings = list(nI1 = nI1, nI0 = nI0, nE = nE, p = p,
                    rho = rho, cov_shift = cov_shift,
                    shape = shape, lambda = lambda,
                    target_cens = target_cens, rate_c = rate_c)
  )
}
