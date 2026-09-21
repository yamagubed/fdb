# Drift curve runner and one-stop study wrapper ==========================

#' Evaluate operating characteristics across a set of drift values
#'
#' For each drift value in \code{drift_set}, runs a simulation under
#' the supplied \code{theta0} (e.g. 0 for type I error, \code{log(0.8)}
#' for power) and the given tuning parameters, then returns a stacked
#' summary across drift values. The reference internal-control sample
#' size used for ESS is taken as \code{nI1 + nI0} from
#' \code{scenario_base}.
#'
#' @param theta0 True treatment effect (log HR) used in the
#'   simulations.
#' @param drift_set Numeric vector of drift values (log HR).
#' @param scenario_base A scenario list. \code{theta0} and \code{delta0}
#'   are overwritten internally for each drift value.
#' @param lambdas Tuning list (as in \code{\link{run_simulation}}).
#' @param nsim Number of replicates per drift value.
#' @param alpha Nominal level.
#' @param seed RNG seed (per-drift offsets are added internally).
#' @param parallel,ncores,robust,eps,n_grid_opt As in
#'   \code{\link{run_simulation}}.
#' @return A data frame stacking the per-drift summaries with an added
#'   \code{driftHR} column.
#' @export
run_drift_curve <- function(theta0,
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
                            n_grid_opt = DEFAULT_N_GRID_OPT) {

  out_list <- vector("list", length(drift_set))

  for (k in seq_along(drift_set)) {
    sc        <- scenario_base
    sc$theta0 <- theta0
    sc$delta0 <- drift_set[k]

    simres <- run_simulation(
      nsim = nsim, scenario = sc, lambdas = lambdas,
      alpha = alpha, one_sided = TRUE,
      seed = seed + 10000L * k,
      parallel = parallel, ncores = ncores,
      robust = robust, eps = eps, n_grid_opt = n_grid_opt
    )

    NS <- sc$nI1 + sc$nI0
    simres <- add_ess_to_simulation_result(simres, NS = NS,
                                           ref_method = "InternalOnly")

    tmp <- simres$summary
    tmp$theta0  <- theta0
    tmp$delta0  <- sc$delta0
    tmp$driftHR <- exp(sc$delta0)
    tmp$nsim    <- nsim

    out_list[[k]] <- tmp
  }

  do.call(rbind, out_list)
}

#' One-stop wrapper: calibrate lambda, then evaluate type I and power
#'
#' Optionally runs design-stage lambda calibration for all five
#' borrowing methods, then evaluates type I error and power curves
#' across the supplied drift set using the calibrated lambdas. If
#' \code{do_calibration = FALSE}, default lambdas are used.
#'
#' @param scenario_base A scenario list. Defaults to
#'   \code{\link{scenario_S1}}.
#' @param drift_hr_range Numeric vector of length 2 giving the
#'   inclusive range of drift HR values for the final curves.
#' @param drift_by_hr Step in HR units for the final drift grid.
#' @param drift_hr_values_cal Numeric vector of drift HR values used
#'   for the (reduced) calibration drift grid. Defaults to a
#'   five-point grid.
#' @param alpha Nominal one-sided level for the final analysis.
#' @param alt_hr Alternative-hypothesis HR used in the power curve
#'   (e.g. 0.8).
#' @param do_calibration Logical; if \code{TRUE}, run lambda calibration
#'   before the final curves.
#' @param lambda_grid Initial coarse lambda grid for calibration.
#' @param lambda_grid_fine_length Number of points in the fine
#'   refinement grid.
#' @param nsim_cal,nsim_confirm Calibration replicates per drift
#'   value.
#' @param nsim_curve Replicates per drift value for the final curves.
#' @param alpha_cal Calibration threshold (defaults to \code{alpha}).
#' @param two_stage_cal Use two-stage calibration.
#' @param confirm_full_drift Run the full-drift confirmation stage.
#' @param primary_calibration_inference Which inference type is used
#'   for tuning during calibration (\code{"sandwich"} or
#'   \code{"model_based"}).
#' @param early_stop_cal,cal_stop_rule,cal_select_rule Drift-level
#'   early stopping and selection rules.
#' @param parallel,ncores Parallelization controls.
#' @param robust,eps,n_grid_opt Estimation controls.
#' @param seed RNG seed.
#' @param rho_mcp MCP transition fraction in (0, 1), default 0.1.
#' @param export_dir If non-\code{NULL}, write CSV outputs and an RDS
#'   of metadata to this directory.
#' @return A list with \code{scenario_base}, \code{drift_hr_range},
#'   \code{drift_set}, \code{drift_set_cal}, \code{alpha},
#'   \code{alpha_cal}, \code{alt_hr}, \code{lambdas} (final tuning),
#'   \code{calibration} (calibration output, or \code{NULL}),
#'   \code{type1_curve}, and \code{power_curve}.
#' @export
run_fdb_study <- function(
  scenario_base    = scenario_S1,
  drift_hr_range   = c(0.8, 1.2),
  drift_by_hr      = 0.05,
  drift_hr_values_cal = c(0.80, 0.90, 1.00, 1.10, 1.20),
  alpha            = 0.025,
  alt_hr           = 0.8,
  do_calibration   = TRUE,
  lambda_grid      = exp(seq(log(0.02), log(2.0), length.out = 6)),
  lambda_grid_fine_length = 6,
  nsim_cal         = 500,
  nsim_confirm     = nsim_cal,
  nsim_curve       = 1000,
  alpha_cal        = alpha,
  two_stage_cal    = TRUE,
  confirm_full_drift = TRUE,
  primary_calibration_inference = c("sandwich", "model_based"),
  early_stop_cal   = TRUE,
  cal_stop_rule    = c("point", "upper95"),
  cal_select_rule  = c("point", "upper95"),
  parallel         = FALSE,
  ncores           = NULL,
  robust           = FALSE,
  eps              = SMOOTH_EPS,
  n_grid_opt       = DEFAULT_N_GRID_OPT,
  seed             = 1,
  export_dir       = NULL,
  rho_mcp          = DEFAULT_RHO_MCP
) {
  cal_stop_rule   <- match.arg(cal_stop_rule)
  cal_select_rule <- match.arg(cal_select_rule)
  primary_calibration_inference <- match.arg(primary_calibration_inference)

  # Full drift grid used for final type I/power curves and, optionally,
  # for the confirmation stage of calibration.
  drift_set <- make_drift_set(drift_hr_range = drift_hr_range,
                              by_hr = drift_by_hr)

  # Reduced drift grid used in the faster calibration stages.
  if (is.null(drift_hr_values_cal)) {
    drift_set_cal <- drift_set
  } else {
    drift_set_cal <- make_drift_set_from_values(drift_hr_values_cal,
                                                include_zero = TRUE)
  }

  lambdas_use <- list(
    lambda_li = 0.20, gamma_li = 1,
    lambda_p1 = 0.20,
    lambda_p2 = 0.20, gate_c = 1.64, gate_tau = 0.25,
    lambda_p3 = 0.20, gamma_mcp = 3, rho_mcp = rho_mcp,
    lambda_p4 = 0.20,
    delta_bounds = DEFAULT_DELTA_BOUNDS
  )

  calib <- NULL
  if (do_calibration) {
    calib <- calibrate_all_lambdas(
      lambda_grid       = lambda_grid,
      scenario_base     = scenario_base,
      drift_set         = drift_set,
      drift_set_cal     = drift_set_cal,
      drift_set_confirm = if (confirm_full_drift) drift_set else drift_set_cal,
      nsim_cal          = nsim_cal,
      nsim_confirm      = nsim_confirm,
      alpha             = alpha,
      alpha_cal         = alpha_cal,
      seed              = seed + 100,
      parallel          = parallel,
      ncores            = ncores,
      robust            = robust,
      eps               = eps,
      gamma_li          = lambdas_use$gamma_li,
      gate_c            = lambdas_use$gate_c,
      gate_tau          = lambdas_use$gate_tau,
      gamma_mcp         = lambdas_use$gamma_mcp,
      rho_mcp           = lambdas_use$rho_mcp,
      delta_bounds      = lambdas_use$delta_bounds,
      n_grid_opt        = n_grid_opt,
      two_stage         = two_stage_cal,
      n_fine            = lambda_grid_fine_length,
      primary_inference = primary_calibration_inference,
      confirm_full_drift = confirm_full_drift,
      early_stop_drift  = early_stop_cal,
      stop_rule         = cal_stop_rule,
      select_rule       = cal_select_rule
    )

    lstar <- calib$lambda_star
    get_lstar <- function(m, inf) {
      v <- lstar[m, inf]
      if (!is.finite(v)) {
        warning("Calibration failed for ", m, "/", inf,
                "; using default lambda = 0.20.")
        0.20
      } else {
        v
      }
    }

    lambdas_use$lambda_li <- get_lstar("Li", primary_calibration_inference)
    lambdas_use$lambda_p1 <- get_lstar("P1", primary_calibration_inference)
    lambdas_use$lambda_p2 <- get_lstar("P2", primary_calibration_inference)
    lambdas_use$lambda_p3 <- get_lstar("P3", primary_calibration_inference)
    lambdas_use$lambda_p4 <- get_lstar("P4", primary_calibration_inference)
  }

  type1_curve <- run_drift_curve(
    theta0 = 0, drift_set = drift_set,
    scenario_base = scenario_base, lambdas = lambdas_use,
    nsim = nsim_curve, alpha = alpha,
    seed = seed + 1000,
    parallel = parallel, ncores = ncores,
    robust = robust, eps = eps, n_grid_opt = n_grid_opt
  )

  power_curve <- run_drift_curve(
    theta0 = log(alt_hr), drift_set = drift_set,
    scenario_base = scenario_base, lambdas = lambdas_use,
    nsim = nsim_curve, alpha = alpha,
    seed = seed + 2000,
    parallel = parallel, ncores = ncores,
    robust = robust, eps = eps, n_grid_opt = n_grid_opt
  )

  if (!is.null(export_dir)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
    utils::write.csv(type1_curve,
                     file.path(export_dir, "type1_curve.csv"),
                     row.names = FALSE)
    utils::write.csv(power_curve,
                     file.path(export_dir, "power_curve.csv"),
                     row.names = FALSE)
    if (!is.null(calib)) {
      utils::write.csv(calib$calibration_table,
                       file.path(export_dir, "lambda_calibration_table.csv"),
                       row.names = FALSE)
      utils::write.csv(calib$calibration_summary,
                       file.path(export_dir, "lambda_calibration_summary.csv"),
                       row.names = FALSE)
      utils::write.csv(as.data.frame(as.table(calib$lambda_star)),
                       file.path(export_dir, "lambda_star.csv"),
                       row.names = FALSE)
    }
    meta <- list(drift_hr_range = drift_hr_range,
                 drift_by_hr = drift_by_hr,
                 drift_set = drift_set,
                 drift_hr_values_cal = drift_hr_values_cal,
                 drift_set_cal = drift_set_cal,
                 alpha = alpha, alpha_cal = alpha_cal, alt_hr = alt_hr,
                 nsim_cal = nsim_cal, nsim_confirm = nsim_confirm,
                 nsim_curve = nsim_curve,
                 lambda_grid = lambda_grid,
                 lambda_grid_fine_length = lambda_grid_fine_length,
                 two_stage_cal = two_stage_cal,
                 confirm_full_drift = confirm_full_drift,
                 primary_calibration_inference = primary_calibration_inference,
                 parallel = parallel, ncores = ncores,
                 robust = robust, eps = eps,
                 n_grid_opt = n_grid_opt, seed = seed,
                 early_stop_cal = early_stop_cal,
                 cal_stop_rule = cal_stop_rule,
                 cal_select_rule = cal_select_rule,
                 scenario_base = scenario_base,
                 lambdas_use = lambdas_use)
    saveRDS(meta, file.path(export_dir, "run_metadata.rds"))
  }

  list(
    scenario_base  = scenario_base,
    drift_hr_range = drift_hr_range,
    drift_set      = drift_set,
    drift_set_cal  = drift_set_cal,
    alpha          = alpha,
    alpha_cal      = alpha_cal,
    alt_hr         = alt_hr,
    lambdas        = lambdas_use,
    calibration    = calib,
    type1_curve    = type1_curve,
    power_curve    = power_curve
  )
}
