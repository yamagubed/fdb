# Smoothed penalty functions and curvatures ==============================
#
# For a smoothed L1 penalty p_lambda,eps(delta) = c_lambda * |delta|_eps,
# where |delta|_eps = sqrt(delta^2 + eps^2), we have
#   d/d delta |delta|_eps = delta / sqrt(delta^2 + eps^2)
#   d^2/d delta^2 |delta|_eps = eps^2 / (delta^2 + eps^2)^(3/2)
# The curvature at the origin is 1/eps (finite).
#
# This file collects the smoothed absolute-value primitive, the L1, MCP,
# and logistic-gate building blocks, and their analytic curvatures used
# for sandwich variance estimation. Stabilization (flooring curvature
# at zero) is applied selectively for non-convex or evidence-gated
# penalties whose raw analytic curvature may be negative in finite
# samples.
# ----------------------------------------------------------------------

#' Smoothed absolute value and its derivatives
#'
#' \code{|delta|_eps = sqrt(delta^2 + eps^2)} and its first two
#' derivatives, used to render L1-type borrowing penalties twice
#' differentiable for sandwich variance estimation.
#'
#' @param delta Scalar drift value.
#' @param eps Smoothing parameter (positive, typically 1e-3).
#' @return Numeric scalar.
#' @keywords internal
#' @noRd
smooth_abs <- function(delta, eps = SMOOTH_EPS) {
  sqrt(delta^2 + eps^2)
}

#' @rdname smooth_abs
#' @keywords internal
#' @noRd
smooth_abs_d1 <- function(delta, eps = SMOOTH_EPS) {
  delta / sqrt(delta^2 + eps^2)
}

#' @rdname smooth_abs
#' @keywords internal
#' @noRd
smooth_abs_d2 <- function(delta, eps = SMOOTH_EPS) {
  eps^2 / (delta^2 + eps^2)^(3 / 2)
}

#' Smoothed L1 penalty
#'
#' @param delta Scalar drift value.
#' @param w Effective penalty weight (\code{c_lambda}).
#' @param eps Smoothing parameter.
#' @return The smoothed L1 penalty value.
#' @keywords internal
#' @noRd
pen_L1 <- function(delta, w, eps = SMOOTH_EPS) {
  w * smooth_abs(delta, eps)
}

#' Smoothed minimax concave penalty
#'
#' @param delta Scalar drift value.
#' @param lambda MCP \code{lambda} (effective penalty strength).
#' @param gamma MCP shape parameter (\code{> 1}).
#' @param eps Smoothing parameter.
#' @return The smoothed MCP value.
#' @keywords internal
#' @noRd
pen_MCP <- function(delta, lambda, gamma = 3, eps = SMOOTH_EPS) {
  ad <- abs(delta)
  if (ad <= gamma * lambda) {
    lambda * smooth_abs(delta, eps) - (delta^2) / (2 * gamma)
  } else {
    0.5 * gamma * lambda^2
  }
}

#' Logistic evidence gate
#'
#' \code{g(t) = (1 + exp((t - c) / tau))^{-1}}.
#'
#' @param t Standardized drift statistic.
#' @param c Evidence threshold.
#' @param tau Smoothness parameter.
#' @return Numeric scalar in (0, 1).
#' @keywords internal
#' @noRd
gate_logistic <- function(t, c = 1.64, tau = 0.25) {
  1 / (1 + exp((t - c) / tau))
}

# Curvatures used in sandwich bread matrix -------------------------------

#' Curvature of a smoothed L1 penalty
#'
#' For \code{p_lambda,eps(delta) = c_lambda * |delta|_eps}, returns
#' \code{c_lambda * eps^2 / (delta^2 + eps^2)^(3/2)}.
#'
#' @keywords internal
#' @noRd
pen_curvature_L1 <- function(delta_hat, c_lambda, eps = SMOOTH_EPS) {
  c_lambda * smooth_abs_d2(delta_hat, eps)
}

#' Curvature of the smooth evidence-gated L1 penalty (full product rule)
#'
#' @keywords internal
#' @noRd
pen_curvature_gated_L1 <- function(delta_hat, lambda, se_delta,
                                   c = 1.64, tau = 0.25,
                                   eps = SMOOTH_EPS) {
  r   <- smooth_abs(delta_hat, eps)
  rp  <- smooth_abs_d1(delta_hat, eps)
  rpp <- smooth_abs_d2(delta_hat, eps)

  se_delta <- pmax(se_delta, NUMERIC_ZERO)
  t   <- r / se_delta
  g   <- gate_logistic(t, c = c, tau = tau)
  gp  <- -(1 / tau) * g * (1 - g)
  gpp <- (1 / tau^2) * g * (1 - g) * (1 - 2 * g)

  f_r  <- g + r * gp / se_delta
  f_rr <- 2 * gp / se_delta + r * gpp / (se_delta^2)
  lambda * (f_rr * rp^2 + f_r * rpp)
}

#' Curvature of the smoothed minimax concave penalty
#'
#' The active region uses
#' \code{lambda_eff * |delta|_eps - delta^2 / (2 * gamma_mcp)}
#' (smoothed near the origin); the tail region returns zero curvature.
#'
#' @keywords internal
#' @noRd
pen_curvature_MCP <- function(delta_hat, lambda_eff, gamma_mcp,
                              eps = SMOOTH_EPS) {
  ad <- abs(delta_hat)
  if (ad <= gamma_mcp * lambda_eff) {
    lambda_eff * smooth_abs_d2(delta_hat, eps) - 1 / gamma_mcp
  } else {
    0
  }
}

#' Stabilize a raw curvature value
#'
#' Used for non-convex or evidence-gated penalties whose analytic
#' curvature can be negative in finite samples. Non-finite values are
#' replaced by zero; if \code{floor_zero = TRUE}, negative values are
#' floored at zero as a finite-sample regularization device.
#'
#' @keywords internal
#' @noRd
stabilize_curvature <- function(curv, floor_zero = TRUE) {
  if (!is.finite(curv)) return(0)
  if (floor_zero) return(max(curv, 0))
  curv
}
