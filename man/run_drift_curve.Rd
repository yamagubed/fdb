\name{run_drift_curve}
\alias{run_drift_curve}
\title{Evaluate operating characteristics across a set of drift values}
\usage{
run_drift_curve(theta0,
                drift_set,
                scenario_base,
                lambdas,
                nsim = 1000,
                alpha = 0.025,
                seed = 1,
                parallel = FALSE,
                ncores = NULL,
                robust = FALSE,
                eps = SMOOTH_EPS,
                n_grid_opt = DEFAULT_N_GRID_OPT)
}
\arguments{
\item{theta0}{True treatment effect (log HR) used in the simulations.}
\item{drift_set}{Numeric vector of drift values on the log-HR scale.}
\item{scenario_base}{A scenario list. \code{theta0} and \code{delta0}
  are overwritten internally for each drift value.}
\item{lambdas}{Tuning list (as in \code{\link{run_simulation}}).}
\item{nsim}{Number of replicates per drift value.}
\item{alpha}{Nominal level.}
\item{seed}{RNG seed (per-drift offsets are added internally).}
\item{parallel, ncores, robust, eps, n_grid_opt}{As in
  \code{\link{run_simulation}}.}
}
\value{
A data frame stacking the per-drift summaries with an added
\code{driftHR} column.
}
\description{
For each drift value, runs a simulation under the supplied
\code{theta0} (e.g. 0 for type I error, \code{log(0.8)} for power)
and the given tuning parameters, then returns a stacked summary
across drift values. The reference internal-control sample size used
for ESS is taken as \code{nI1 + nI0} from \code{scenario_base}.
}
