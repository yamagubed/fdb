\name{calibrate_all_lambdas}
\alias{calibrate_all_lambdas}
\title{Calibrate lambda for all borrowing methods}
\usage{
calibrate_all_lambdas(lambda_grid,
                      scenario_base,
                      drift_set,
                      drift_set_cal = NULL,
                      drift_set_confirm = drift_set,
                      nsim_cal = 500,
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
                      two_stage = TRUE,
                      n_fine = 6,
                      primary_inference = c("sandwich", "model_based"),
                      confirm_full_drift = TRUE,
                      early_stop_drift = TRUE,
                      stop_rule = c("point", "upper95"),
                      select_rule = c("point", "upper95"))
}
\arguments{
\item{lambda_grid}{Initial (coarse) lambda grid used for all methods.}
\item{scenario_base}{A scenario list.}
\item{drift_set}{The full drift set on which to evaluate the final
  confirmation (also used as the default calibration drift set).}
\item{drift_set_cal}{Reduced drift set for the calibration stages.}
\item{drift_set_confirm}{Drift set for the optional confirmation
  stage (defaults to \code{drift_set}).}
\item{nsim_cal, nsim_confirm}{Replicates per drift value in
  calibration and confirmation stages.}
\item{alpha, alpha_cal, seed, parallel, ncores, robust, eps, gamma_li, gate_c, gate_tau, gamma_mcp, delta_bounds, n_grid_opt, n_fine, primary_inference, confirm_full_drift, early_stop_drift, stop_rule, select_rule}{
  As in \code{\link{calibrate_lambda_grid_two_stage}}.}
\item{two_stage}{Logical; use two-stage calibration if \code{TRUE},
  otherwise single-stage.}
}
\value{
A list with \code{calibration_table} (long form),
\code{calibration_summary} (per-method, per-lambda worst case),
\code{lambda_star} (5 x 2 matrix indexed by method and inference),
\code{lambda_star_table} (long form), and \code{method_calibrations}
(per-method full output).
}
\description{
Runs lambda calibration for all five penalized borrowing methods
(Li adaptive lasso plus P1-P4) and returns a unified set of
calibrated lambdas under both model-based and sandwich inference.
}
