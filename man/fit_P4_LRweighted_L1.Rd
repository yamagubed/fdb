\name{fit_P4_LRweighted_L1}
\alias{fit_P4_LRweighted_L1}
\title{Likelihood-ratio-weighted L1 penalty (P4)}
\usage{
fit_P4_LRweighted_L1(dat, xnames = NULL, lambda,
                     delta_bounds = DEFAULT_DELTA_BOUNDS,
                     full_fit = NULL, nodelta_fit = NULL,
                     robust = FALSE, eps = SMOOTH_EPS)
}
\arguments{
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions.}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected.}
\item{lambda}{Penalty strength.}
\item{delta_bounds}{Optimization interval for \eqn{\delta}.}
\item{full_fit}{Optional pre-computed full Cox fit.}
\item{nodelta_fit}{Optional pre-computed restricted (no-Z) Cox fit.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
\item{eps}{Smoothing parameter for \eqn{|\delta|_\varepsilon}.}
}
\value{
A list of estimates and inference quantities, including the
likelihood ratio statistic \code{LR0} and effective weight \code{w}.
}
\description{
Implements
\eqn{p_\lambda(\delta) = \lambda |\delta| \exp\{-\tfrac{1}{2} LR(\delta = 0)\}},
where \eqn{LR(\delta = 0)} is the likelihood ratio statistic for
testing the null of no population drift. Borrowing is downweighted
when external data conflict strongly with internal controls.
}
