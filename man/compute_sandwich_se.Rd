\name{compute_sandwich_se}
\alias{compute_sandwich_se}
\title{Sandwich standard error for a penalized Cox borrowing estimator}
\usage{
compute_sandwich_se(dat, delta_hat, theta_hat, beta_hat,
                    pen_curv, xnames)
}
\arguments{
\item{dat}{A data frame with columns \code{time}, \code{status},
  \code{T}, \code{Z}, and covariates named in \code{xnames}.}
\item{delta_hat}{Penalized estimate of the drift parameter.}
\item{theta_hat}{Penalized estimate of the treatment effect.}
\item{beta_hat}{Penalized estimate of the covariate coefficient vector.}
\item{pen_curv}{Penalty curvature
  \eqn{p''_{\lambda,\varepsilon}(\hat\delta)}.}
\item{xnames}{Character vector of covariate names.}
}
\value{
A list with components \code{se_sand}, \code{Var_sand},
\code{coef_names}, \code{theta_idx}, \code{delta_idx},
\code{pen_curv}, and \code{note}.
}
\description{
Computes the plug-in sandwich variance for the smoothed penalized
Cox M-estimator. The Cox information and score residuals are
evaluated at the penalized parameter estimate, and the penalty
curvature is added to the (delta, delta) element of the bread.
}
\details{
The penalized estimator satisfies
\deqn{\Psi_\varepsilon(\eta) = U(\eta) - e_\delta\, p'_{\lambda,\varepsilon}(\delta) = 0,}
where \eqn{U(\eta)} is the Cox partial-likelihood score and
\eqn{e_\delta} selects the drift coordinate. A local plug-in sandwich
approximation is
\eqn{\mathrm{Var}(\hat\eta) = A_\varepsilon^{-1} \hat B A_\varepsilon^{-T}}
with \eqn{A_\varepsilon = I(\hat\eta) + p''_{\lambda,\varepsilon}(\hat\delta) e_\delta e_\delta^T}
and \eqn{\hat B = \sum_i \hat v_i(\hat\eta) \hat v_i(\hat\eta)^T},
where \eqn{\hat v_i} denotes the subject-level Cox score residual.

This is a local plug-in approximation holding first-stage estimates fixed.
It does not establish conditional or unconditional variance validity.
The meat uses subject-level Cox score residuals evaluated at the penalized
coefficients. Clipped curvature modifies only the variance calculation.
}
