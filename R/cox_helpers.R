# Cox model fitting helpers ==============================================

#' Fit a profiled Cox model with delta as a fixed offset
#'
#' Fits a Cox proportional hazards model on \code{(T, X)} with
#' \code{offset(delta * Z)}. Used during profiling optimization
#' over \code{delta}.
#'
#' @param dat A data frame with columns \code{time}, \code{status},
#'   \code{T}, \code{Z}, and covariates named in \code{xnames}.
#' @param delta Scalar drift value to use as the offset multiplier.
#' @param xnames Character vector of covariate names.
#' @param robust If \code{TRUE}, use robust (Lin-Wei) variance from
#'   \code{survival::coxph}.
#' @return A list with \code{theta_hat}, \code{se_theta},
#'   \code{beta_hat}, \code{loglik}, and the fitted \code{coxph} object.
#' @keywords internal
#' @noRd
cox_fit_profile <- function(dat, delta, xnames, robust = FALSE) {
  fml <- stats::as.formula(paste0("survival::Surv(time, status) ~ T + ",
                                  paste(xnames, collapse = " + "),
                                  " + offset(off)"))
  dat$off <- delta * dat$Z
  fit <- survival::coxph(fml, data = dat, ties = "efron", robust = robust)
  s <- summary(fit)
  se_col <- if (robust && "robust se" %in% colnames(s$coef)) {
    "robust se"
  } else {
    "se(coef)"
  }
  list(
    theta_hat = as.numeric(stats::coef(fit)["T"]),
    se_theta  = as.numeric(s$coef["T", se_col]),
    beta_hat  = as.numeric(stats::coef(fit)[xnames]),
    loglik    = fit$loglik[2],
    fit       = fit
  )
}

#' Fit the unpenalized full Cox model on (T, Z, X)
#'
#' Returns the unpenalized estimator of the drift parameter together
#' with the treatment effect, used by data-adaptive penalties as the
#' first-stage estimator.
#'
#' @inheritParams cox_fit_profile
#' @return A list of estimates and the fitted \code{coxph} object.
#' @keywords internal
#' @noRd
cox_fit_full <- function(dat, xnames, robust = FALSE) {
  fml <- stats::as.formula(paste0("survival::Surv(time, status) ~ T + Z + ",
                                  paste(xnames, collapse = " + ")))
  fit <- survival::coxph(fml, data = dat, ties = "efron", robust = robust)
  s <- summary(fit)
  se_col <- if (robust && "robust se" %in% colnames(s$coef)) {
    "robust se"
  } else {
    "se(coef)"
  }
  list(
    theta_hat = as.numeric(stats::coef(fit)["T"]),
    se_theta  = as.numeric(s$coef["T", se_col]),
    delta_hat = as.numeric(stats::coef(fit)["Z"]),
    se_delta  = as.numeric(s$coef["Z", se_col]),
    loglik    = fit$loglik[2],
    fit       = fit
  )
}

#' Fit a Cox model omitting the drift term Z
#'
#' Used as the restricted estimator (under \eqn{\delta = 0}) for the
#' likelihood-ratio-weighted penalty and for full-pooling reference
#' analyses.
#'
#' @inheritParams cox_fit_profile
#' @return A list of estimates and the fitted \code{coxph} object.
#' @keywords internal
#' @noRd
cox_fit_nodelta <- function(dat, xnames, robust = FALSE) {
  fml <- stats::as.formula(paste0("survival::Surv(time, status) ~ T + ",
                                  paste(xnames, collapse = " + ")))
  fit <- survival::coxph(fml, data = dat, ties = "efron", robust = robust)
  s <- summary(fit)
  se_col <- if (robust && "robust se" %in% colnames(s$coef)) {
    "robust se"
  } else {
    "se(coef)"
  }
  list(
    theta_hat = as.numeric(stats::coef(fit)["T"]),
    se_theta  = as.numeric(s$coef["T", se_col]),
    loglik    = fit$loglik[2],
    fit       = fit
  )
}
