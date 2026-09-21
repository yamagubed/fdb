\name{calibrate_lambda_grid}
\alias{calibrate_lambda_grid}
\title{Single-stage lambda calibration for one borrowing method}
\usage{
calibrate_lambda_grid(method = c("Li", "P1", "P2", "P3", "P4"),
                      lambda_grid,
                      scenario_base,
                      drift_set,
                      nsim = 300,
                      alpha = 0.025,
                      alpha_cal = alpha,
                      seed = 1,
                      parallel = FALSE,
                      ncores = NULL,
                      robust = FALSE,
                      eps = SMOOTH_EPS,
                      gamma_li = 1,
                      gate_c = 1.64,
                      gate_tau = 0.25,
                      gamma_mcp = 3,
                      delta_bounds = DEFAULT_DELTA_BOUNDS,
                      n_grid_opt = DEFAULT_N_GRID_OPT,
                      early_stop_drift = FALSE,
                      stop_rule = c("point", "upper95"),
                      select_rule = c("point", "upper95"),
  rho_mcp = DEFAULT_RHO_MCP)
}
\arguments{
\item{rho_mcp}{MCP transition fraction in (0, 1). The draft software
  default is 0.1 and is not a calibrated manuscript value. The transition
  half-width is rho_mcp * gamma_mcp * lambda_eff.}
\item{method}{One of \code{"Li"}, \code{"P1"}, \code{"P2"},
  \code{"P3"}, \code{"P4"}.}
\item{lambda_grid}{Numeric vector of candidate lambda values.}
\item{scenario_base}{A scenario list as accepted by
  \code{\link{run_simulation}}. Its \code{theta0} and \code{delta0}
  are overwritten internally.}
\item{drift_set}{Numeric vector of drift values (log HR) on which to
  evaluate type I error.}
\item{nsim}{Number of replicates per drift value.}
\item{alpha}{Nominal level used for rejection decisions.}
\item{alpha_cal}{Calibration threshold (worst-case rejection rate
  must not exceed this).}
\item{seed}{RNG seed.}
\item{parallel}{Logical; enable parallel evaluation across
  replicates.}
\item{ncores}{Number of cores for parallel evaluation.}
\item{robust}{Use robust (Lin-Wei) Cox SEs.}
\item{eps}{Smoothing parameter.}
\item{gamma_li}{Adaptive lasso exponent.}
\item{gate_c, gate_tau}{P2 gate parameters.}
\item{gamma_mcp}{MCP shape parameter for P3.}
\item{delta_bounds}{Optimization interval for \eqn{\delta}.}
\item{n_grid_opt}{Coarse-grid points for non-convex objectives.}
\item{early_stop_drift}{Logical; if \code{TRUE}, drop lambda values
  that have already failed the calibration constraint after
  evaluating a subset of drift values.}
\item{stop_rule}{Either \code{"point"} (use point estimate of
  worst-case rejection rate for early stopping) or \code{"upper95"}
  (use the Monte Carlo 95\% upper confidence bound).}
\item{select_rule}{Selection rule for the final calibrated lambda
  (same options as \code{stop_rule}).}
}
\value{
A list with elements \code{method}, \code{details} (per-lambda x
drift), \code{summary} (per-lambda worst-case),
\code{calibration_table} (long format), and \code{lambda_star}
(selected lambdas for both inference types).
}
\description{
For each lambda in \code{lambda_grid}, simulates \code{nsim}
replicates under each drift value in \code{drift_set} (with
\eqn{\theta_0 = 0}), records the model-based and sandwich rejection
rates, and selects the largest lambda whose worst-case rejection
rate over the drift set does not exceed \code{alpha_cal}.
}
\examples{
\donttest{
cal <- calibrate_lambda_grid(
  method = "P1",
  lambda_grid = exp(seq(log(0.05), log(1), length.out = 4)),
  scenario_base = scenario_S1,
  drift_set = make_drift_set_from_values(c(0.9, 1.0, 1.1)),
  nsim = 100, seed = 1)
cal$lambda_star
}
}
