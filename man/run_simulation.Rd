\name{run_simulation}
\alias{run_simulation}
\title{Run a Monte Carlo simulation under a fixed scenario}
\usage{
run_simulation(nsim = 200,
               scenario,
               lambdas,
               alpha = 0.025,
               one_sided = TRUE,
               seed = 1,
               parallel = FALSE,
               ncores = NULL,
               robust = FALSE,
               eps = SMOOTH_EPS,
               n_grid_opt = DEFAULT_N_GRID_OPT)
}
\arguments{
\item{nsim}{Number of Monte Carlo replicates.}
\item{scenario}{A list of data-generating parameters with elements
  \code{nI1}, \code{nI0}, \code{nE}, \code{theta0}, \code{delta0},
  \code{p}, \code{beta}, \code{rho}, \code{cov_shift}, \code{shape},
  \code{lambda}, and \code{target_cens}. See \code{\link{scenario_S1}}
  for an example.}
\item{lambdas}{A list of tuning parameters with elements
  \code{lambda_li}, \code{gamma_li}, \code{lambda_p1},
  \code{lambda_p2}, \code{gate_c}, \code{gate_tau},
  \code{lambda_p3}, \code{gamma_mcp}, \code{lambda_p4}, and
  \code{delta_bounds}. See \code{\link{lambdas_default}} for an example.}
\item{alpha}{Nominal one-sided (or two-sided) significance level.}
\item{one_sided}{Logical; if \code{TRUE}, rejection corresponds to a
  hazard ratio less than one.}
\item{seed}{RNG seed.}
\item{parallel}{Logical; if \code{TRUE}, parallelize replicates using
  \code{parallel::parLapply}.}
\item{ncores}{Number of cores for parallel execution.}
\item{robust}{Use robust (Lin-Wei) Cox standard errors.}
\item{eps}{Smoothing parameter for \eqn{|\delta|_\varepsilon}.}
\item{n_grid_opt}{Number of grid points for the coarse search in
  non-convex objectives.}
}
\value{
A list with components \code{raw} (per-replicate results),
\code{summary} (aggregated for both inference types, with an
\code{inference} column), \code{scenario}, \code{lambdas}, and
\code{settings}.
}
\description{
For each of \code{nsim} replicates, simulates a hybrid-control Cox
dataset and fits all borrowing methods. Returns per-replicate raw
results and aggregated summaries under both model-based and
sandwich inference.
}
\examples{
\donttest{
sim_out <- run_simulation(nsim = 50, scenario = scenario_S1,
                          lambdas = lambdas_default, alpha = 0.025, seed = 1)
subset(sim_out$summary, inference == "sandwich")
}
}
