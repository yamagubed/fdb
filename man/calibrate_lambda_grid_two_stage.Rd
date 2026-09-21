\name{calibrate_lambda_grid_two_stage}
\alias{calibrate_lambda_grid_two_stage}
\title{Two-stage lambda calibration for one borrowing method}
\usage{
calibrate_lambda_grid_two_stage(method = c("Li", "P1", "P2", "P3", "P4"),
                                lambda_grid_coarse,
                                scenario_base,
                                drift_set_cal,
                                drift_set_confirm = NULL,
                                nsim_cal = 300,
                                nsim_confirm = nsim_cal,
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
                                n_fine = 6,
                                primary_inference = c("sandwich", "model_based"),
                                confirm_full_drift = TRUE,
                                early_stop_drift = TRUE,
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
\item{lambda_grid_coarse}{Initial coarse lambda grid.}
\item{scenario_base}{A scenario list.}
\item{drift_set_cal}{Reduced drift set for calibration.}
\item{drift_set_confirm}{Deprecated for final lambda selection; retained
  for backward compatibility.}
\item{nsim_cal}{Replicates per drift value in the calibration stages.}
\item{nsim_confirm}{Replicates per drift value in the confirmation
  stage.}
\item{alpha, alpha_cal}{Nominal level and calibration threshold.}
\item{seed}{RNG seed (subsequent stages use offsets of this seed).}
\item{parallel, ncores}{Parallelization controls.}
\item{robust, eps, gamma_li, gate_c, gate_tau, gamma_mcp, delta_bounds, n_grid_opt}{
  Same as in \code{\link{calibrate_lambda_grid}}.}
\item{n_fine}{Number of refined points, before adding the original coarse grid.}
\item{primary_inference}{Which inference type drives the refinement
  (default \code{"sandwich"}).}
\item{confirm_full_drift}{Logical; perform the confirmation stage on
  \code{drift_set_cal}.}
\item{early_stop_drift, stop_rule, select_rule}{As in
  \code{\link{calibrate_lambda_grid}}.}
}
\value{
A list with components \code{method}, \code{coarse}, \code{fine},
\code{confirm}, \code{calibration_summary}, \code{calibration_table},
\code{lambda_star}, and \code{final_stage}.
}
\description{
Performs a coarse-grid calibration on a reduced drift set
(\code{drift_set_cal}), refines the grid around the calibrated lambda,
repeats the fine-grid calibration on \code{drift_set_cal}, and
optionally confirms candidates on the same drift set, \code{drift_set_cal}.
The fine-stage grid is the union of the original coarse grid and the
refined points. Smaller candidates are retained for confirmation, even
if they failed an earlier stage. Selection uses the final stage's results;
a passing coarse-stage result is never substituted for confirmation.
}
