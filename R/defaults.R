# Default scenarios and tuning parameters ================================

#' Example simulation scenario
#'
#' A reference scenario used in the examples and unit tests. Provides
#' 150 internal treated, 150 internal concurrent control, and 300
#' external control subjects, with 5 correlated covariates and a
#' Weibull baseline.
#'
#' @format A list with elements \code{nI1}, \code{nI0}, \code{nE},
#'   \code{theta0}, \code{delta0}, \code{p}, \code{beta}, \code{rho},
#'   \code{cov_shift}, \code{shape}, \code{lambda}, and
#'   \code{target_cens}.
#' @export
scenario_S1 <- list(
  nI1 = 150, nI0 = 150, nE = 300,
  theta0 = 0, delta0 = 0,
  p = 5, beta = rep(0.2, 5),
  rho = 0.3,
  cov_shift = rep(0, 5),
  shape = 1.2, lambda = 0.02,
  target_cens = 0.2
)

#' Default tuning parameters
#'
#' A reference tuning list for the borrowing methods, used when no
#' calibrated tuning is supplied. The default penalty strengths
#' (\eqn{\lambda = 0.20}) are reasonable starting values for the
#' default scenario but should generally be replaced by calibrated
#' values for confirmatory use.
#'
#' @format A list with elements \code{lambda_li}, \code{gamma_li},
#'   \code{lambda_p1}, \code{lambda_p2}, \code{gate_c}, \code{gate_tau},
#'   \code{lambda_p3}, \code{gamma_mcp}, \code{rho_mcp}, \code{lambda_p4}, and
#'   \code{delta_bounds}.
#' @export
lambdas_default <- list(
  lambda_li = 0.20, gamma_li = 1,
  lambda_p1 = 0.20,
  lambda_p2 = 0.20, gate_c = 1.64, gate_tau = 0.25,
  lambda_p3 = 0.20, gamma_mcp = 3, rho_mcp = 0.1,
  lambda_p4 = 0.20,
  delta_bounds = c(-2, 2)
)
