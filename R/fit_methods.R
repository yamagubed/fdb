# Method-specific fit functions ==========================================
#
# Each method returns a list with at least the components
#   method, theta_hat, se_theta (model-based), se_sand (sandwich),
#   z, z_sand, delta_hat, pen_curv,
# together with method-specific quantities (effective penalty weight,
# evidence-gate value, likelihood ratio statistic, etc.).
# ----------------------------------------------------------------------

#' Internal-control-only Cox analysis
#'
#' Fits a Cox proportional hazards model on the internal (Z = 0)
#' subjects only, without any borrowing from external controls.
#'
#' @param dat A data frame from \code{\link{simulate_hybrid_cox}} or
#'   conforming to its column conventions.
#' @param xnames Character vector of covariate names. If \code{NULL},
#'   automatically detected as columns matching \code{"^X\\d+$"}.
#' @param robust Use robust (Lin-Wei) Cox standard errors.
#' @return A list of estimates and inference quantities.
#' @export
fit_internal_only <- function(dat, xnames = NULL, robust = FALSE) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  fit <- cox_fit_nodelta(dat[dat$Z == 0, , drop = FALSE], xnames,
                         robust = robust)
  z <- fit$theta_hat / fit$se_theta
  list(
    method    = "InternalOnly",
    theta_hat = fit$theta_hat,
    se_theta  = fit$se_theta,
    se_sand   = fit$se_theta,
    z         = z,
    z_sand    = z,
    delta_hat = NA_real_,
    pen_curv  = 0
  )
}

#' Naive pooled Cox analysis
#'
#' Fits a Cox proportional hazards model treating concurrent and
#' external controls as exchangeable (full pooling), without any
#' indicator for external-control membership.
#'
#' @inheritParams fit_internal_only
#' @return A list of estimates and inference quantities.
#' @export
fit_naive_pooled <- function(dat, xnames = NULL, robust = FALSE) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  fit <- cox_fit_nodelta(dat, xnames, robust = robust)
  z <- fit$theta_hat / fit$se_theta
  list(
    method    = "NaivePooled",
    theta_hat = fit$theta_hat,
    se_theta  = fit$se_theta,
    se_sand   = fit$se_theta,
    z         = z,
    z_sand    = z,
    delta_hat = NA_real_,
    pen_curv  = 0
  )
}

#' Adaptive lasso borrowing of Li et al. (2023)
#'
#' Implements the adaptive lasso borrowing penalty
#' \eqn{p_\lambda(\delta) = \lambda |\hat\delta_0|^{-\gamma} |\delta|},
#' where \eqn{\hat\delta_0} is the unpenalized maximum partial
#' likelihood estimator of the drift parameter.
#'
#' @param dat A data frame from \code{\link{simulate_hybrid_cox}} or
#'   conforming to its column conventions.
#' @param xnames Character vector of covariate names. If \code{NULL},
#'   automatically detected.
#' @param lambda Penalty strength (\eqn{\lambda \ge 0}).
#' @param gamma Adaptive-weight exponent (typically 1).
#' @param delta_bounds Numeric vector of length 2 giving the
#'   optimization interval for \eqn{\delta}.
#' @param full_fit Optional pre-computed full Cox fit (output of an
#'   internal helper). If \code{NULL}, the fit is computed.
#' @param robust Use robust (Lin-Wei) Cox standard errors.
#' @param eps Smoothing parameter for \eqn{|\delta|_\varepsilon}.
#' @return A list of estimates and inference quantities, including
#'   the effective penalty weight \code{w}, the initial estimator
#'   \code{delta0_hat}, and its standard error \code{se_delta}.
#'
#' @references
#' Li, R., Lin, R., Huang, J., Tian, L., and Zhu, J. (2023).
#' A frequentist approach to dynamic borrowing.
#' \emph{Biometrical Journal} 65(7), 2100406.
#'
#' @export
fit_li_adaptive_lasso <- function(dat, xnames = NULL, lambda, gamma = 1,
                                  delta_bounds = DEFAULT_DELTA_BOUNDS,
                                  full_fit = NULL, robust = FALSE,
                                  eps = SMOOTH_EPS) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  if (is.null(full_fit)) {
    full_fit <- cox_fit_full(dat, xnames, robust = robust)
  }
  delta0_hat <- full_fit$delta_hat
  if (abs(delta0_hat) < NUMERIC_ZERO) {
    warning("delta0_hat near zero; adaptive lasso weight may be unstable")
  }

  w <- lambda * pmax(abs(delta0_hat), NUMERIC_ZERO)^(-gamma)

  obj <- function(delta) {
    prof <- cox_fit_profile(dat, delta, xnames, robust = robust)
    (-prof$loglik) + pen_L1(delta, w, eps = eps)
  }

  sol <- optimize_delta(dat, xnames, obj, delta_bounds, robust = robust)

  pen_curv <- pen_curvature_L1(sol$delta_hat, c_lambda = w, eps = eps)
  pen_curv <- stabilize_curvature(pen_curv)

  build_method_result("LiAdaptiveLasso", sol, dat, xnames, pen_curv,
                      extra = list(w = w,
                                   delta0_hat = delta0_hat,
                                   se_delta   = full_fit$se_delta))
}

#' Precision-weighted L1 penalty (P1)
#'
#' Implements
#' \eqn{p_\lambda(\delta) = \lambda |\delta| / \widehat{SE}(\hat\delta_0)},
#' where \eqn{\widehat{SE}(\hat\delta_0)} is the standard error of the
#' unpenalized drift estimator.
#'
#' @inheritParams fit_li_adaptive_lasso
#' @return A list of estimates and inference quantities, including the
#'   effective penalty weight \code{w} and the first-stage standard
#'   error \code{se_delta}.
#' @export
fit_P1_precision_L1 <- function(dat, xnames = NULL, lambda,
                                delta_bounds = DEFAULT_DELTA_BOUNDS,
                                full_fit = NULL, robust = FALSE,
                                eps = SMOOTH_EPS) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  if (is.null(full_fit)) {
    full_fit <- cox_fit_full(dat, xnames, robust = robust)
  }
  w <- lambda / pmax(full_fit$se_delta, NUMERIC_ZERO)

  obj <- function(delta) {
    prof <- cox_fit_profile(dat, delta, xnames, robust = robust)
    (-prof$loglik) + pen_L1(delta, w, eps = eps)
  }

  sol <- optimize_delta(dat, xnames, obj, delta_bounds, robust = robust)
  pen_curv <- pen_curvature_L1(sol$delta_hat, c_lambda = w, eps = eps)
  pen_curv <- stabilize_curvature(pen_curv)

  build_method_result("P1_SEScaledL1", sol, dat, xnames, pen_curv,
                      extra = list(se_delta = full_fit$se_delta, w = w))
}

#' Smoothed integrated-gate penalty (P2)
#'
#' Implements
#' \deqn{p_{\lambda,\varepsilon}(\delta) = \lambda\int_\varepsilon^{\sqrt{\delta^2+\varepsilon^2}} [1+\exp\{(u/\widehat{SE}(\hat\delta_0)-c)/\tau\}]^{-1}du,}
#' which provides a continuous relaxation of test-then-pool procedures
#' driven by standardized evidence against commensurability.
#'
#' @inheritParams fit_li_adaptive_lasso
#' @param c Evidence threshold (default 1.64, one-sided 5\% gate).
#' @param tau Gate smoothness parameter (positive).
#' @param n_grid_opt Number of grid points for the coarse search in
#'   the non-convex objective.
#' @return A list of estimates and inference quantities, including the
#'   evidence-gate value \code{g_hat} and the raw (unstabilized)
#'   penalty curvature \code{pen_curv_raw}.
#' @export
fit_P2_gated_L1 <- function(dat, xnames = NULL, lambda,
                            c = 1.64, tau = 0.25,
                            delta_bounds = DEFAULT_DELTA_BOUNDS,
                            full_fit = NULL, robust = FALSE,
                            eps = SMOOTH_EPS,
                            n_grid_opt = DEFAULT_N_GRID_OPT) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  if (is.null(full_fit)) {
    full_fit <- cox_fit_full(dat, xnames, robust = robust)
  }
  se_d <- full_fit$se_delta
  validate_gate(lambda, se_d, c, tau, eps)

  obj <- function(delta) {
    prof <- cox_fit_profile(dat, delta, xnames, robust = robust)
    (-prof$loglik) + pen_integrated_gate(delta, lambda, se_d, c, tau, eps)
  }

  sol <- optimize_delta(dat, xnames, obj, delta_bounds,
                       robust = robust, n_grid = n_grid_opt)

  t_hat <- smooth_abs(sol$delta_hat, eps) / se_d
  g_hat <- gate_logistic(t_hat, c = c, tau = tau)
  pen_curv_raw <- pen_curvature_gated_L1(sol$delta_hat, lambda = lambda,
                                         se_delta = se_d,
                                         c = c, tau = tau, eps = eps)
  pen_curv <- stabilize_curvature(pen_curv_raw)

  build_method_result("P2_GatedL1", sol, dat, xnames, pen_curv,
                      extra = list(se_delta = full_fit$se_delta,
                                   lambda = lambda, c = c, tau = tau,
                                   g_hat = g_hat,
                                   pen_curv_raw = pen_curv_raw))
}

#' Information-adaptive minimax concave penalty (P3)
#'
#' Implements
#' \deqn{p_\lambda(\delta) = MCP(\delta;\,\lambda/\widehat{SE}(\hat\delta_0),\,\gamma_{MCP}),}
#' where MCP is the minimax concave penalty of Zhang (2010). Reduces
#' bias when moderate population drift is present while retaining
#' strong shrinkage near \eqn{\delta = 0}.
#'
#' @inheritParams fit_li_adaptive_lasso
#' @param gamma_mcp MCP shape parameter (\code{> 1}).
#' @param rho_mcp MCP transition fraction in (0, 1); h = rho_mcp *
#'   gamma_mcp * lambda_eff. The software default 0.1 is not calibrated.
#' @details MCP is smoothed at both the origin and the flat-tail transition
#'   by integrating the continuously differentiable slope described in the
#'   manuscript. First-stage quantities, including h, are held fixed.
#' @param n_grid_opt Number of grid points for the coarse search in
#'   the non-convex objective.
#' @return A list of estimates and inference quantities, including the
#'   effective MCP \code{lambda_eff} and the raw curvature
#'   \code{pen_curv_raw}.
#'
#' @references
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under
#' minimax concave penalty. \emph{The Annals of Statistics} 38(2),
#' 894-942.
#'
#' @export
fit_P3_info_MCP <- function(dat, xnames = NULL, lambda, gamma_mcp = 3,
                            delta_bounds = DEFAULT_DELTA_BOUNDS,
                            full_fit = NULL, robust = FALSE,
                            eps = SMOOTH_EPS,
                            n_grid_opt = DEFAULT_N_GRID_OPT,
                            rho_mcp = DEFAULT_RHO_MCP) {
  validate_mcp(lambda, gamma_mcp, eps, rho_mcp)
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  if (is.null(full_fit)) {
    full_fit <- cox_fit_full(dat, xnames, robust = robust)
  }
  validate_penalty_scalar(full_fit$se_delta, "se_delta", 0, TRUE)
  lam_eff <- lambda / full_fit$se_delta

  obj <- function(delta) {
    prof <- cox_fit_profile(dat, delta, xnames, robust = robust)
    (-prof$loglik) + pen_MCP(delta, lambda = lam_eff, gamma = gamma_mcp,
                             eps = eps, rho = rho_mcp)
  }

  sol <- optimize_delta(dat, xnames, obj, delta_bounds,
                       robust = robust, n_grid = n_grid_opt)

  pen_curv_raw <- pen_curvature_MCP(sol$delta_hat, lambda_eff = lam_eff,
                                    gamma_mcp = gamma_mcp, eps = eps,
                                    rho = rho_mcp)
  pen_curv <- stabilize_curvature(pen_curv_raw)

  build_method_result("P3_SEScaledMCP", sol, dat, xnames, pen_curv,
                      extra = list(se_delta = full_fit$se_delta,
                                   lambda_eff = lam_eff,
                                   gamma_mcp = gamma_mcp,
                                   rho_mcp = rho_mcp,
                                   h = rho_mcp * gamma_mcp * lam_eff,
                                   pen_curv_raw = pen_curv_raw))
}

#' Likelihood-ratio-weighted L1 penalty (P4)
#'
#' Implements
#' \eqn{p_\lambda(\delta) = \lambda |\delta| \exp\{-\tfrac{1}{2} LR(\delta = 0)\}},
#' where \eqn{LR(\delta = 0)} is the likelihood ratio statistic for
#' testing the null of no population drift. Borrowing is downweighted
#' when external data conflict strongly with internal controls.
#'
#' @inheritParams fit_li_adaptive_lasso
#' @param nodelta_fit Optional pre-computed restricted (no-Z) Cox fit.
#' @return A list of estimates and inference quantities, including the
#'   likelihood ratio statistic \code{LR0} and effective weight \code{w}.
#' @export
fit_P4_LRweighted_L1 <- function(dat, xnames = NULL, lambda,
                                 delta_bounds = DEFAULT_DELTA_BOUNDS,
                                 full_fit = NULL, nodelta_fit = NULL,
                                 robust = FALSE, eps = SMOOTH_EPS) {
  if (is.null(xnames)) {
    xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  }
  if (is.null(full_fit))    full_fit    <- cox_fit_full(dat, xnames, robust = robust)
  if (is.null(nodelta_fit)) nodelta_fit <- cox_fit_nodelta(dat, xnames, robust = robust)

  LR0 <- 2 * (full_fit$loglik - nodelta_fit$loglik)
  w   <- lambda * exp(-0.5 * LR0)

  obj <- function(delta) {
    prof <- cox_fit_profile(dat, delta, xnames, robust = robust)
    (-prof$loglik) + pen_L1(delta, w, eps = eps)
  }

  sol      <- optimize_delta(dat, xnames, obj, delta_bounds, robust = robust)
  pen_curv <- pen_curvature_L1(sol$delta_hat, c_lambda = w, eps = eps)
  pen_curv <- stabilize_curvature(pen_curv)

  build_method_result("P4_LRWeightedL1", sol, dat, xnames, pen_curv,
                      extra = list(LR0 = LR0, w = w))
}
