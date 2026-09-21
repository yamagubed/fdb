\name{fit_P3_info_MCP}
\alias{fit_P3_info_MCP}
\title{Information-adaptive minimax concave penalty (P3)}
\usage{
fit_P3_info_MCP(dat, xnames = NULL, lambda, gamma_mcp = 3,
                delta_bounds = DEFAULT_DELTA_BOUNDS,
                full_fit = NULL, robust = FALSE,
                eps = SMOOTH_EPS, n_grid_opt = DEFAULT_N_GRID_OPT,
  rho_mcp = DEFAULT_RHO_MCP)
}
\arguments{
\item{rho_mcp}{MCP transition fraction in (0, 1). The draft software
  default is 0.1 and is not a calibrated manuscript value. The transition
  half-width is rho_mcp * gamma_mcp * lambda_eff.}
\item{dat}{A data frame from \code{\link{simulate_hybrid_cox}} or
  conforming to its column conventions.}
\item{xnames}{Character vector of covariate names. If \code{NULL},
  automatically detected.}
\item{lambda}{Penalty strength.}
\item{gamma_mcp}{MCP shape parameter (\code{> 1}).}
\item{delta_bounds}{Optimization interval for \eqn{\delta}.}
\item{full_fit}{Optional pre-computed full Cox fit.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
\item{eps}{Smoothing parameter for \eqn{|\delta|_\varepsilon}.}
\item{n_grid_opt}{Number of grid points for the coarse search in the
  non-convex objective.}
}
\value{
A list of estimates and inference quantities, including the
effective MCP \code{lambda_eff} and the raw curvature
\code{pen_curv_raw}.
}
\description{
Implements
\deqn{p_\lambda(\delta) = MCP(\delta;\,\lambda/\widehat{SE}(\hat\delta_0),\,\gamma_{MCP}),}
where MCP is the minimax concave penalty of Zhang (2010). Reduces
bias when moderate population drift is present while retaining strong
shrinkage near \eqn{\delta = 0}.
}
\references{
Zhang, C.-H. (2010). Nearly unbiased variable selection under
minimax concave penalty. \emph{The Annals of Statistics} 38(2),
894-942.
}

\details{
The smooth penalty integrates the manuscript C1 MCP slope from zero to
\eqn{\sqrt{\delta^2+\varepsilon^2}-\varepsilon}. Both the origin and
the flat-tail transition are smoothed. The returned \code{rho_mcp} and
\code{h} record the transition settings. Raw curvature is retained;
negative curvature is clipped only for the variance calculation.
}
