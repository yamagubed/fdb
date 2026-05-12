\name{fit_li_adaptive_lasso}
\alias{fit_li_adaptive_lasso}
\title{Adaptive lasso borrowing of Li et al. (2023)}
\usage{
fit_li_adaptive_lasso(dat, xnames = NULL, lambda, gamma = 1,
                      delta_bounds = DEFAULT_DELTA_BOUNDS,
                      full_fit = NULL, robust = FALSE,
                      eps = SMOOTH_EPS)
}
\arguments{
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions.}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected.}
\item{lambda}{Penalty strength (\eqn{\lambda \ge 0}).}
\item{gamma}{Adaptive-weight exponent (typically 1).}
\item{delta_bounds}{Numeric vector of length 2 giving the
  optimization interval for \eqn{\delta}.}
\item{full_fit}{Optional pre-computed full Cox fit. If \code{NULL},
  the fit is computed internally.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
\item{eps}{Smoothing parameter for \eqn{|\delta|_\varepsilon}.}
}
\value{
A list of estimates and inference quantities, including the
effective penalty weight \code{w}, the initial estimator
\code{delta0_hat}, and its standard error \code{se_delta}.
}
\description{
Implements the adaptive lasso borrowing penalty
\eqn{p_\lambda(\delta) = \lambda |\hat\delta_0|^{-\gamma} |\delta|},
where \eqn{\hat\delta_0} is the unpenalized maximum partial
likelihood estimator of the drift parameter.
}
\references{
Li, R., Lin, R., Huang, J., Tian, L., and Zhu, J. (2023).
A frequentist approach to dynamic borrowing.
\emph{Biometrical Journal} 65(7), 2100406.
}
