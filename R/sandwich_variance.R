# Sandwich variance for the penalized Cox M-estimator ===================
#
# Under the penalized likelihood formulation, the treatment effect
# estimator theta_hat is obtained jointly with the drift estimator
# delta_hat and nuisance covariate effects beta_hat. The penalized
# estimator satisfies
#
#   Psi_eps(eta) = U(eta) - e_delta * p'_lambda,eps(delta) = 0
#
# where U(eta) is the Cox partial-likelihood score, e_delta selects the
# delta coordinate, and p_lambda,eps(.) is the smoothed penalty.
#
# A local plug-in sandwich approximation is
#
#   Var(eta_hat) = A_eps^{-1} * B_hat * A_eps^{-T}
#
# with A_eps = I(eta_hat) + p''_lambda,eps(delta_hat) e_delta e_delta^T
# and B_hat = sum_i vhat_i(eta_hat) vhat_i(eta_hat)^T. Below we use the total-
# scale equivalent of this expression to avoid duplicate (1/n) factors.
#
# Quantities used to define penalty weights or gates (delta_hat_0,
# SE(delta_hat_0), LR statistic) are treated as fixed after their
# first-stage estimation. This approximation is
# not guaranteed valid for either conditional or unconditional variance.
# The vhat_i are subject-level Cox score residuals, not event contributions.
# Clipping changes the variance calculation only, not the fitted coefficients.
# ----------------------------------------------------------------------

#' Sandwich standard error for a penalized Cox borrowing estimator
#'
#' Computes the plug-in sandwich variance for the smoothed penalized
#' Cox M-estimator. The Cox information and score residuals are
#' evaluated at the penalized parameter estimate, and the penalty
#' curvature is added to the (delta, delta) element of the bread.
#' First-stage estimates are held fixed; this does not establish conditional
#' or unconditional variance validity. Clipped curvature modifies the
#' variance calculation without changing the fitted coefficients.
#'
#' @param dat A data frame with columns \code{time}, \code{status},
#'   \code{T}, \code{Z}, and covariates named in \code{xnames}.
#' @param delta_hat Penalized estimate of the drift parameter.
#' @param theta_hat Penalized estimate of the treatment effect.
#' @param beta_hat Penalized estimate of the covariate coefficient
#'   vector.
#' @param pen_curv Penalty curvature \eqn{p''_{\lambda,\varepsilon}(\hat\delta)}.
#' @param xnames Character vector of covariate names.
#' @return A list with components:
#'   \describe{
#'     \item{\code{se_sand}}{Sandwich standard error for theta_hat.}
#'     \item{\code{Var_sand}}{Full sandwich variance matrix.}
#'     \item{\code{coef_names}, \code{theta_idx}, \code{delta_idx}}{
#'       Coefficient bookkeeping.}
#'     \item{\code{pen_curv}}{Penalty curvature used in the bread.}
#'     \item{\code{note}}{Status note ("OK" on success).}
#'   }
#' @export
compute_sandwich_se <- function(dat, delta_hat, theta_hat, beta_hat,
                                pen_curv, xnames) {
  fml_full <- stats::reformulate(c("T", "Z", xnames),
                                 response = "survival::Surv(time, status)")

  eta_hat <- c(theta_hat, delta_hat, beta_hat)
  names(eta_hat) <- c("T", "Z", xnames)

  # Evaluate the Cox machinery at eta_hat without iterating.
  eval_fit <- tryCatch(
    survival::coxph(fml_full, data = dat, ties = "efron", init = eta_hat,
                    control = survival::coxph.control(iter.max = 0, eps = 1e-9),
                    robust = FALSE, x = TRUE, model = TRUE),
    error = function(e) NULL
  )
  if (is.null(eval_fit)) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "Cox evaluation at penalized eta failed"))
  }

  score_resid <- tryCatch(
    stats::residuals(eval_fit, type = "score"),
    error = function(e) NULL
  )
  if (is.null(score_resid)) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "score residuals unavailable at penalized eta"))
  }
  if (!is.matrix(score_resid)) score_resid <- as.matrix(score_resid)

  # Total-scale meat B_total = sum_i vhat_i(eta_hat) vhat_i(eta_hat)^T.
  B_total <- crossprod(score_resid)

  vcov_eval <- tryCatch(stats::vcov(eval_fit), error = function(e) NULL)
  if (is.null(vcov_eval)) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "vcov unavailable"))
  }

  I_total <- tryCatch(solve(vcov_eval), error = function(e) NULL)
  if (is.null(I_total)) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "information matrix inversion failed"))
  }

  coef_names <- names(stats::coef(eval_fit))
  delta_idx  <- which(coef_names == "Z")
  theta_idx  <- which(coef_names == "T")

  if (length(delta_idx) == 0) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "Z (delta) not found"))
  }
  if (length(theta_idx) == 0) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "T (theta) not found"))
  }

  q <- length(coef_names)
  e_delta <- numeric(q)
  e_delta[delta_idx] <- 1

  # Total-scale bread: I_total + p''_{lambda,eps}(delta_hat) e_delta e_delta^T.
  A_eps <- I_total + pen_curv * outer(e_delta, e_delta)

  A_inv <- tryCatch(solve(A_eps), error = function(e) NULL)
  if (is.null(A_inv)) {
    return(list(se_sand = NA_real_, Var_sand = NULL,
                note = "A_eps inversion failed"))
  }

  Var_sand <- A_inv %*% B_total %*% t(A_inv)
  se_sand  <- sqrt(max(0, Var_sand[theta_idx, theta_idx]))

  list(
    se_sand    = se_sand,
    Var_sand   = Var_sand,
    coef_names = coef_names,
    theta_idx  = theta_idx,
    delta_idx  = delta_idx,
    pen_curv   = pen_curv,
    note       = "OK"
  )
}
