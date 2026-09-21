# Higher-level fit dispatchers ===========================================

#' Fit a single penalized borrowing method (cached version)
#'
#' Internal helper: fits one penalized method given pre-computed
#' full and (optionally) no-delta Cox fits. Used in calibration where
#' the first-stage fits are reused across many lambda values.
#'
#' @keywords internal
#' @noRd
fit_one_penalized_method_cached <- function(dat,
                                            xnames,
                                            full_fit,
                                            nodelta_fit = NULL,
                                            method = c("Li", "P1", "P2", "P3", "P4"),
                                            lambda,
                                            gamma_li = 1,
                                            gate_c = 1.64,
                                            gate_tau = 0.25,
                                            gamma_mcp = 3,
                                            delta_bounds = DEFAULT_DELTA_BOUNDS,
                                            robust = FALSE,
                                            eps = SMOOTH_EPS,
                                            n_grid_opt = DEFAULT_N_GRID_OPT,
                                            rho_mcp = DEFAULT_RHO_MCP) {
  method <- match.arg(method)

  if (method == "Li") {
    return(fit_li_adaptive_lasso(
      dat, xnames,
      lambda = lambda,
      gamma = gamma_li,
      delta_bounds = delta_bounds,
      full_fit = full_fit,
      robust = robust,
      eps = eps
    ))
  }

  if (method == "P1") {
    return(fit_P1_precision_L1(
      dat, xnames,
      lambda = lambda,
      delta_bounds = delta_bounds,
      full_fit = full_fit,
      robust = robust,
      eps = eps
    ))
  }

  if (method == "P2") {
    return(fit_P2_gated_L1(
      dat, xnames,
      lambda = lambda,
      c = gate_c,
      tau = gate_tau,
      delta_bounds = delta_bounds,
      full_fit = full_fit,
      robust = robust,
      eps = eps,
      n_grid_opt = n_grid_opt
    ))
  }

  if (method == "P3") {
    return(fit_P3_info_MCP(
      dat, xnames,
      lambda = lambda,
      gamma_mcp = gamma_mcp, rho_mcp = rho_mcp,
      delta_bounds = delta_bounds,
      full_fit = full_fit,
      robust = robust,
      eps = eps,
      n_grid_opt = n_grid_opt
    ))
  }

  if (method == "P4") {
    if (is.null(nodelta_fit)) {
      nodelta_fit <- cox_fit_nodelta(dat, xnames, robust = robust)
    }
    return(fit_P4_LRweighted_L1(
      dat, xnames,
      lambda = lambda,
      delta_bounds = delta_bounds,
      full_fit = full_fit,
      nodelta_fit = nodelta_fit,
      robust = robust,
      eps = eps
    ))
  }
}

#' Fit a single penalized borrowing method
#'
#' Dispatches to the appropriate user-facing fit function for one of
#' the implemented penalized borrowing methods.
#'
#' @param dat A data frame conforming to the column conventions of
#'   \code{\link{simulate_hybrid_cox}}.
#' @param method One of \code{"Li"}, \code{"P1"}, \code{"P2"},
#'   \code{"P3"}, \code{"P4"}.
#' @param lambda Penalty strength.
#' @param gamma_li Exponent for the adaptive lasso weight (Li method).
#' @param gate_c,gate_tau Threshold and smoothness for the P2 gate.
#' @param gamma_mcp MCP shape parameter for P3.
#' @param rho_mcp MCP transition fraction in (0, 1), default 0.1.
#' @param delta_bounds Optimization interval for \code{delta}.
#' @param robust Use robust (Lin-Wei) Cox standard errors.
#' @param eps Smoothing parameter for \eqn{|\delta|_\varepsilon}.
#' @param n_grid_opt Number of grid points for the coarse search
#'   (used by P2 and P3).
#' @return A list of estimates and inference quantities (same format
#'   as the underlying fit functions).
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_cox(nI1 = 100, nI0 = 100, nE = 200,
#'                            theta0 = log(0.8), delta0 = 0)
#' fit_one_penalized_method(sim$data, method = "P1", lambda = 0.2)
#' @export
fit_one_penalized_method <- function(dat,
                                     method = c("Li", "P1", "P2", "P3", "P4"),
                                     lambda,
                                     gamma_li = 1,
                                     gate_c = 1.64,
                                     gate_tau = 0.25,
                                     gamma_mcp = 3,
                                     delta_bounds = DEFAULT_DELTA_BOUNDS,
                                     robust = FALSE,
                                     eps = SMOOTH_EPS,
                                     n_grid_opt = DEFAULT_N_GRID_OPT,
                                     rho_mcp = DEFAULT_RHO_MCP) {
  method <- match.arg(method)
  xnames <- grep("^X\\d+$", names(dat), value = TRUE)
  full_fit <- cox_fit_full(dat, xnames, robust = robust)
  nodelta_fit <- if (method == "P4") {
    cox_fit_nodelta(dat, xnames, robust = robust)
  } else {
    NULL
  }

  fit_one_penalized_method_cached(
    dat = dat, xnames = xnames, full_fit = full_fit,
    nodelta_fit = nodelta_fit,
    method = method, lambda = lambda, gamma_li = gamma_li,
    gate_c = gate_c, gate_tau = gate_tau, gamma_mcp = gamma_mcp, rho_mcp = rho_mcp,
    delta_bounds = delta_bounds, robust = robust, eps = eps,
    n_grid_opt = n_grid_opt
  )
}

#' Fit all borrowing methods on a single dataset
#'
#' Convenience wrapper that fits internal-only, naive pooled, Li
#' adaptive lasso, and all four likelihood-informed penalties (P1-P4)
#' on a single hybrid-control dataset and returns a tidy data frame of
#' results. Both model-based and sandwich standard errors are
#' reported.
#'
#' @param dat A data frame conforming to the column conventions of
#'   \code{\link{simulate_hybrid_cox}}.
#' @param lambda_li,gamma_li Tuning for the adaptive lasso method.
#' @param lambda_p1 Tuning for the precision-weighted L1 (P1).
#' @param lambda_p2,gate_c,gate_tau Tuning for the smoothed integrated-gate
#'   L1 (P2).
#' @param rho_mcp MCP transition fraction in (0, 1), default 0.1.
#' @param lambda_p3,gamma_mcp Tuning for the information-adaptive MCP (P3).
#' @param lambda_p4 Tuning for the likelihood-ratio-weighted L1 (P4).
#' @param delta_bounds Optimization interval for \code{delta}.
#' @param robust Use robust (Lin-Wei) Cox standard errors.
#' @param eps Smoothing parameter for \eqn{|\delta|_\varepsilon}.
#' @param n_grid_opt Number of grid points for the coarse search
#'   (used by P2 and P3).
#' @return A data frame with one row per method, columns
#'   \code{method}, \code{theta_hat}, \code{se_theta} (model-based),
#'   \code{se_sand} (sandwich), \code{z}, \code{z_sand},
#'   \code{delta_hat}, \code{pen_curv}, \code{pen_curv_raw}.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_cox(nI1 = 100, nI0 = 100, nE = 200,
#'                            theta0 = log(0.8), delta0 = 0)
#' fit_all_methods(sim$data)
#' @export
fit_all_methods <- function(dat,
                            lambda_li = 0.2, gamma_li = 1,
                            lambda_p1 = 0.2,
                            lambda_p2 = 0.2, gate_c = 1.64, gate_tau = 0.25,
                            lambda_p3 = 0.2, gamma_mcp = 3,
                            lambda_p4 = 0.2,
                            delta_bounds = DEFAULT_DELTA_BOUNDS,
                            robust = FALSE,
                            eps = SMOOTH_EPS,
                            n_grid_opt = DEFAULT_N_GRID_OPT,
                            rho_mcp = DEFAULT_RHO_MCP) {
  xnames <- grep("^X\\d+$", names(dat), value = TRUE)

  # Cache expensive unpenalized fits
  full_fit    <- cox_fit_full(dat, xnames, robust = robust)
  nodelta_fit <- cox_fit_nodelta(dat, xnames, robust = robust)

  out <- list()
  out[[1]] <- fit_internal_only(dat, xnames, robust = robust)
  out[[2]] <- fit_naive_pooled(dat, xnames, robust = robust)
  out[[3]] <- fit_li_adaptive_lasso(dat, xnames, lambda = lambda_li,
                                    gamma = gamma_li,
                                    delta_bounds = delta_bounds,
                                    full_fit = full_fit,
                                    robust = robust, eps = eps)
  out[[4]] <- fit_P1_precision_L1(dat, xnames, lambda = lambda_p1,
                                  delta_bounds = delta_bounds,
                                  full_fit = full_fit,
                                  robust = robust, eps = eps)
  out[[5]] <- fit_P2_gated_L1(dat, xnames, lambda = lambda_p2,
                              c = gate_c, tau = gate_tau,
                              delta_bounds = delta_bounds,
                              full_fit = full_fit,
                              robust = robust, eps = eps,
                              n_grid_opt = n_grid_opt)
  out[[6]] <- fit_P3_info_MCP(dat, xnames, lambda = lambda_p3,
                              gamma_mcp = gamma_mcp, rho_mcp = rho_mcp,
                              delta_bounds = delta_bounds,
                              full_fit = full_fit,
                              robust = robust, eps = eps,
                              n_grid_opt = n_grid_opt)
  out[[7]] <- fit_P4_LRweighted_L1(dat, xnames, lambda = lambda_p4,
                                   delta_bounds = delta_bounds,
                                   full_fit = full_fit,
                                   nodelta_fit = nodelta_fit,
                                   robust = robust, eps = eps)

  res <- do.call(rbind, lapply(out, function(x) {
    data.frame(
      method       = x$method,
      theta_hat    = x$theta_hat,
      se_theta     = x$se_theta,
      se_sand      = if (!is.null(x$se_sand))  x$se_sand  else NA_real_,
      z            = if (!is.null(x$z))        x$z        else x$theta_hat / x$se_theta,
      z_sand       = if (!is.null(x$z_sand))   x$z_sand   else NA_real_,
      delta_hat    = if (!is.null(x$delta_hat)) x$delta_hat else NA_real_,
      pen_curv     = if (!is.null(x$pen_curv))  x$pen_curv  else NA_real_,
      pen_curv_raw = if (!is.null(x$pen_curv_raw)) {
        x$pen_curv_raw
      } else if (!is.null(x$pen_curv)) {
        x$pen_curv
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  res
}
