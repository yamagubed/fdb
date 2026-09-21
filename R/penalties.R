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

# Validate mathematical domains before profile optimization can catch errors.
validate_penalty_scalar <- function(x, name, lower = 0, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    stop(name, if (strict) " must be finite and > " else " must be finite and >= ", lower)
  }
}

validate_mcp <- function(lambda, gamma, eps, rho) {
  validate_penalty_scalar(lambda, "lambda")
  validate_penalty_scalar(gamma, "gamma_mcp", 1, TRUE)
  validate_penalty_scalar(eps, "eps", 0, TRUE)
  validate_penalty_scalar(rho, "rho_mcp", 0, TRUE)
  if (rho >= 1) stop("rho_mcp must be < 1")
}

validate_gate <- function(lambda, se_delta, c, tau, eps) {
  validate_penalty_scalar(lambda, "lambda")
  validate_penalty_scalar(se_delta, "se_delta", 0, TRUE)
  validate_penalty_scalar(c, "c", 0, TRUE)
  validate_penalty_scalar(tau, "tau", 0, TRUE)
  validate_penalty_scalar(eps, "eps", 0, TRUE)
}

# Logistic gate and integrated P2 penalty. Stable closed-form integration
# avoids numerical quadrature inside each profile-likelihood evaluation.
gate_logistic <- function(t, c = 1.64, tau = 0.25) {
  stats::plogis((c - t) / tau)
}

pen_integrated_gate <- function(delta, lambda, se_delta, c = 1.64,
                                tau = 0.25, eps = SMOOTH_EPS) {
  if (lambda == 0) return(0)
  r <- smooth_abs(delta, eps)
  distance <- delta^2 / (r + eps) # r - eps, without cancellation near zero
  width <- se_delta * tau
  d <- distance / width
  x0 <- (c - eps / se_delta) / tau
  x1 <- x0 - d
  integral <- if (d < 0.1) {
    -log1p(stats::plogis(x0) * expm1(-d))
  } else if (x1 > 0) {
    d + log1p(exp(-x0)) - log1p(exp(-x1))
  } else {
    max(x0, 0) + log1p(exp(-abs(x0))) - log1p(exp(x1))
  }
  lambda * width * integral
}

pen_derivative_gated_L1 <- function(delta, lambda, se_delta,
                                   c = 1.64, tau = 0.25, eps = SMOOTH_EPS) {
  lambda * gate_logistic(smooth_abs(delta, eps) / se_delta, c, tau) *
    smooth_abs_d1(delta, eps)
}

pen_curvature_L1 <- function(delta_hat, c_lambda, eps = SMOOTH_EPS) {
  c_lambda * smooth_abs_d2(delta_hat, eps)
}

# Retain the historical internal name for compatibility; P2 now integrates
# the gate rather than multiplying the gate by the drift magnitude.
pen_curvature_gated_L1 <- function(delta_hat, lambda, se_delta,
                                   c = 1.64, tau = 0.25, eps = SMOOTH_EPS) {
  r <- smooth_abs(delta_hat, eps)
  g <- gate_logistic(r / se_delta, c, tau)
  lambda * (g * smooth_abs_d2(delta_hat, eps) -
              g * (1 - g) / (se_delta * tau) * (delta_hat / r)^2)
}

# Integral, slope and slope derivative of the C1 MCP transition function.
# The returned integral is exactly flat beyond b+h and zero at u=0.
mcp_transition <- function(u, lambda, gamma, rho) {
  if (lambda == 0) return(list(value = 0, slope = 0, derivative = 0))
  b <- gamma * lambda
  h <- rho * b
  left <- b - h
  if (u <= left) {
    return(list(value = lambda * u - u^2 / (2 * gamma),
                slope = lambda - u / gamma, derivative = -1 / gamma))
  }
  if (u < b + h) {
    x <- u - left
    return(list(value = lambda * left - left^2 / (2 * gamma) +
                  h * x / gamma - x^2 / (2 * gamma) + x^3 / (12 * gamma * h),
                slope = (b + h - u)^2 / (4 * gamma * h),
                derivative = -(b + h - u) / (2 * gamma * h)))
  }
  list(value = gamma * lambda^2 / 2 + h^2 / (6 * gamma),
       slope = 0, derivative = 0)
}

pen_MCP <- function(delta, lambda, gamma = 3, eps = SMOOTH_EPS,
                    rho = DEFAULT_RHO_MCP) {
  r <- smooth_abs(delta, eps)
  mcp_transition(delta^2 / (r + eps), lambda, gamma, rho)$value
}

pen_derivative_MCP <- function(delta, lambda, gamma = 3, eps = SMOOTH_EPS,
                               rho = DEFAULT_RHO_MCP) {
  r <- smooth_abs(delta, eps)
  mcp_transition(delta^2 / (r + eps), lambda, gamma, rho)$slope * delta / r
}

pen_curvature_MCP <- function(delta_hat, lambda_eff, gamma_mcp,
                              eps = SMOOTH_EPS, rho = DEFAULT_RHO_MCP) {
  r <- smooth_abs(delta_hat, eps)
  q <- mcp_transition(delta_hat^2 / (r + eps), lambda_eff, gamma_mcp, rho)
  q$derivative * (delta_hat / r)^2 + q$slope * smooth_abs_d2(delta_hat, eps)
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
