# Fixed, population-based censoring calibration. No RNG is used here.
.validate_censor_rate <- function(rate) {
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) || rate < 0) {
    stop("censor_rate must be one finite nonnegative number", call. = FALSE)
  }
  invisible(rate)
}

#' Calculate a fixed exponential censoring rate
#'
#' Solves for the rate giving the target marginal censoring probability among
#' internal subjects, weighted by their treatment allocation. Uses numerical
#' integration of the Weibull event-time distribution; normally distributed
#' covariate effects are integrated by 64-point Gaussian quadrature. No random
#' draws or realized trial data are used. The rate does not depend on external
#' drift or external sample size. Calculate separate rates under the null and
#' alternative if the same internal censoring target is desired under both.
#'
#' @param nI1,nI0 Internal treated and control sample sizes.
#' @param theta0 True treatment log hazard ratio.
#' @param p Number of normally distributed baseline covariates.
#' @param beta Covariate coefficients; NULL uses rep(0.2, p).
#' @param rho Equicorrelation of internal covariates.
#' @param shape,lambda Baseline Weibull cumulative hazard is lambda * t^shape.
#' @param target_cens Target internal censoring proportion in [0, 1).
#' @return A nonnegative exponential rate. Zero means no censoring.
#' @examples
#' calibrate_censor_rate(300, 150, p = 0, target_cens = 0.4)
#' @export
calibrate_censor_rate <- function(nI1 = 150, nI0 = 150, theta0 = 0,
                                  p = 5, beta = NULL, rho = 0,
                                  shape = 1.2, lambda = 0.02,
                                  target_cens = 0.2) {
  .validate_count(nI1, "nI1"); .validate_count(nI0, "nI0")
  .validate_count(p, "p", 0L)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(theta0) || !scalar(shape) || shape <= 0 ||
      !scalar(lambda) || lambda <= 0 || !scalar(rho) || abs(rho) > 1 ||
      !scalar(target_cens) || target_cens < 0 || target_cens >= 1) {
    stop("Invalid population parameters or target_cens", call. = FALSE)
  }
  if (p > 1 && (rho <= -1 / (p - 1) || rho >= 1)) {
    stop("rho must give a positive definite covariate correlation matrix")
  }
  if (is.null(beta)) beta <- rep(0.2, p)
  if (!is.numeric(beta) || length(beta) != p || any(!is.finite(beta))) {
    stop("beta must be finite and have length p")
  }
  if (target_cens == 0) return(0)
  variance <- (1 - rho) * sum(beta^2) + rho * sum(beta)^2
  if (variance > 0) {
    # Eigenvalues and squared first eigenvector components integrate N(0,1).
    J <- matrix(0, 64, 64)
    J[cbind(1:63, 2:64)] <- sqrt(1:63)
    J <- J + t(J)
    eig <- eigen(J, symmetric = TRUE)
    nodes <- sqrt(variance) * eig$values
    weights <- eig$vectors[1, ]^2
  } else {
    nodes <- 0; weights <- 1
  }
  allocation <- nI1 / (nI1 + nI0)
  lp <- c(nodes + theta0, nodes)
  weights <- c(allocation * weights, (1 - allocation) * weights)
  # Solve in dimensionless units for numerical stability across baseline scales.
  time_scale <- exp(-lp / shape)
  if (any(!is.finite(time_scale))) stop("Population parameters exceed numerical range")
  prob <- function(r) {
    if (r == 0) return(0)
    stats::integrate(function(x) {
      vapply(x, function(u) {
        sum(weights * (-expm1(-r * u^(1 / shape) * time_scale))) * exp(-u)
      }, numeric(1))
    }, 0, Inf, rel.tol = 1e-9, abs.tol = 1e-10, subdivisions = 500L)$value
  }
  upper <- 1
  while (prob(upper) < target_cens) {
    upper <- upper * 2
    if (!is.finite(upper)) stop("Could not bracket censoring rate")
  }
  root <- stats::uniroot(function(r) prob(r) - target_cens,
                         c(0, upper), tol = 1e-10)$root
  rate <- root * lambda^(1 / shape)
  .validate_censor_rate(rate)
  rate
}

.resolve_censoring <- function(scenario) {
  if (is.null(scenario$censor_rate)) {
    args <- scenario[intersect(names(scenario), names(formals(calibrate_censor_rate)))]
    scenario$censor_rate <- do.call(calibrate_censor_rate, args)
  }
  .validate_censor_rate(scenario$censor_rate)
  scenario
}
