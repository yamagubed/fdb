\name{fit_P2_gated_L1}
\alias{fit_P2_gated_L1}
\title{Smooth evidence-gated L1 penalty (P2)}
\usage{
fit_P2_gated_L1(dat, xnames = NULL, lambda,
                c = 1.64, tau = 0.25,
                delta_bounds = DEFAULT_DELTA_BOUNDS,
                full_fit = NULL, robust = FALSE,
                eps = SMOOTH_EPS, n_grid_opt = DEFAULT_N_GRID_OPT)
}
\arguments{
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions.}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected.}
\item{lambda}{Penalty strength.}
\item{c}{Evidence threshold (default 1.64, one-sided 5\% gate).}
\item{tau}{Gate smoothness parameter (positive).}
\item{delta_bounds}{Optimization interval for \eqn{\delta}.}
\item{full_fit}{Optional pre-computed full Cox fit.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
\item{eps}{Smoothing parameter for \eqn{|\delta|_\varepsilon}.}
\item{n_grid_opt}{Number of grid points for the coarse search in the
  non-convex objective.}
}
\value{
A list of estimates and inference quantities, including the
evidence-gate value \code{g_hat} and the raw (unstabilized) penalty
curvature \code{pen_curv_raw}.
}
\description{
Implements
\deqn{p_\lambda(\delta) = \lambda |\delta| \big[ 1 + \exp((|\delta|/\widehat{SE}(\hat\delta_0) - c)/\tau) \big]^{-1},}
which provides a continuous relaxation of test-then-pool procedures
driven by standardized evidence against commensurability.
}
