# Lambda calibration logic ===============================================
#
# Speed improvements over a naive implementation:
#   * Calibration fits only the target borrowing method per call.
#   * Model-based and sandwich rejection rates are computed from the
#     same simulated trials.
#   * Two-stage calibration: coarse grid, then fine grid near the
#     borrowing boundary.
#   * Optional confirmation on a wider drift set after fine-stage
#     calibration on a reduced drift set.
#   * Within each simulated trial, the unpenalized full Cox fit and
#     no-delta Cox fit are computed once and reused across all lambda
#     values evaluated for that method.
#   * Optional drift-level early stopping removes lambda values that
#     have already exceeded the calibration threshold.
# ----------------------------------------------------------------------

#' Drift set utilities
#'
#' \code{make_drift_set} builds a regular grid of drift values from a
#' range and step in HR units, optionally including the null drift
#' (\code{HR = 1}). \code{make_drift_set_from_values} converts a
#' user-supplied set of drift HR values to log-HR.
#'
#' @param drift_hr_range Numeric vector of length 2 giving the
#'   inclusive range of drift HR values (e.g. \code{c(0.8, 1.2)}).
#' @param by_hr Step size in HR units.
#' @param include_zero If \code{TRUE}, ensure the null drift
#'   (\eqn{\delta = 0}) is included.
#' @return A numeric vector of drift values on the log-HR scale.
#'
#' @name drift_set
#' @export
make_drift_set <- function(drift_hr_range = c(0.8, 1.2),
                           by_hr = 0.05,
                           include_zero = TRUE) {
  stopifnot(length(drift_hr_range) == 2,
            drift_hr_range[1] > 0, drift_hr_range[2] > 0)
  hr_grid <- seq(drift_hr_range[1], drift_hr_range[2], by = by_hr)
  d <- log(hr_grid)
  if (include_zero) d <- sort(unique(c(0, d)))
  d
}

#' @rdname drift_set
#' @param drift_hr_values Numeric vector of drift HR values.
#' @export
make_drift_set_from_values <- function(drift_hr_values,
                                       include_zero = TRUE) {
  stopifnot(is.numeric(drift_hr_values), all(drift_hr_values > 0))
  d <- log(drift_hr_values)
  if (include_zero) d <- sort(unique(c(0, d)))
  d
}

# Internal helpers ------------------------------------------------------

#' Monte Carlo upper confidence bound for a calibration rejection rate
#'
#' @keywords internal
#' @noRd
.calibration_mc_upper <- function(p_hat, nsim, level = 0.95) {
  z <- stats::qnorm((1 + level) / 2)
  se <- sqrt(pmax(p_hat, 0) * pmax(1 - p_hat, 0) / nsim)
  p_hat + z * se
}

#' Select the largest lambda satisfying the calibration constraint
#'
#' @keywords internal
#' @noRd
.select_lambda_star <- function(summary_tbl,
                                alpha_cal,
                                inference = c("model_based", "sandwich"),
                                select_rule = c("point", "upper95")) {
  inference <- match.arg(inference)
  select_rule <- match.arg(select_rule)
  if (nrow(summary_tbl) == 0) return(NA_real_)

  if (inference == "model_based") {
    v <- if (select_rule == "upper95") {
      summary_tbl$worst_type1_model_based_upper95
    } else {
      summary_tbl$worst_type1_model_based
    }
  } else {
    v <- if (select_rule == "upper95") {
      summary_tbl$worst_type1_sandwich_upper95
    } else {
      summary_tbl$worst_type1_sandwich
    }
  }
  ok <- which(is.finite(v) & v <= alpha_cal)
  if (length(ok) == 0) NA_real_ else summary_tbl$lambda[max(ok)]
}

#' Build a refined log-spaced lambda grid around the calibrated lambda
#'
#' @keywords internal
#' @noRd
.make_refined_lambda_grid <- function(lambda_grid,
                                      lambda_star,
                                      n_fine = 6,
                                      expand_if_at_upper = FALSE) {
  lambda_grid <- sort(unique(lambda_grid))
  if (length(lambda_grid) == 1) return(lambda_grid)
  n_fine <- max(2, as.integer(n_fine))

  if (!is.finite(lambda_star)) {
    lower <- lambda_grid[1]
    upper <- lambda_grid[2]
  } else {
    idx <- max(which(lambda_grid <= lambda_star * (1 + 1e-12)))
    if (idx >= length(lambda_grid)) {
      lower <- lambda_grid[max(1, idx - 1)]
      upper <- lambda_grid[idx]
      if (isTRUE(expand_if_at_upper)) {
        upper <- lambda_grid[idx] * (lambda_grid[idx] / lower)
      }
    } else {
      lower <- if (idx <= 1) lambda_grid[1] else lambda_grid[idx - 1]
      upper <- lambda_grid[idx + 1]
    }
  }

  if (!is.finite(lower) || !is.finite(upper) ||
      lower <= 0 || upper <= 0 || lower == upper) {
    return(lambda_grid)
  }
  sort(unique(exp(seq(log(lower), log(upper), length.out = n_fine))))
}

#' Compute model-based and sandwich rejection across a lambda grid for one trial
#'
#' @keywords internal
#' @noRd
.calibration_rejects_for_dataset <- function(dat,
                                             method,
                                             lambda_grid,
                                             zcrit,
                                             gamma_li,
                                             gate_c,
                                             gate_tau,
                                             gamma_mcp,
                                             delta_bounds,
                                             robust,
                                             eps,
                                             n_grid_opt) {
  xnames <- grep("^X\\d+$", names(dat), value = TRUE)

  full_fit <- tryCatch(cox_fit_full(dat, xnames, robust = robust),
                       error = function(e) NULL)
  if (is.null(full_fit)) {
    return(data.frame(lambda = lambda_grid,
                      reject_mod = NA_real_,
                      reject_sand = NA_real_))
  }

  nodelta_fit <- NULL
  if (method == "P4") {
    nodelta_fit <- tryCatch(cox_fit_nodelta(dat, xnames, robust = robust),
                            error = function(e) NULL)
    if (is.null(nodelta_fit)) {
      return(data.frame(lambda = lambda_grid,
                        reject_mod = NA_real_,
                        reject_sand = NA_real_))
    }
  }

  rows <- lapply(lambda_grid, function(lam) {
    fit <- tryCatch(
      fit_one_penalized_method_cached(
        dat = dat, xnames = xnames, full_fit = full_fit,
        nodelta_fit = nodelta_fit,
        method = method, lambda = lam, gamma_li = gamma_li,
        gate_c = gate_c, gate_tau = gate_tau, gamma_mcp = gamma_mcp,
        delta_bounds = delta_bounds, robust = robust, eps = eps,
        n_grid_opt = n_grid_opt
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(lambda = lam,
                        reject_mod = NA_real_,
                        reject_sand = NA_real_))
    }
    data.frame(
      lambda      = lam,
      reject_mod  = as.numeric(is.finite(fit$z) && fit$z < zcrit),
      reject_sand = as.numeric(is.finite(fit$z_sand) && fit$z_sand < zcrit)
    )
  })

  do.call(rbind, rows)
}

#' Evaluate calibration across a lambda grid at one drift value
#'
#' @keywords internal
#' @noRd
.evaluate_calibration_drift_cached <- function(method,
                                               lambda_grid,
                                               scenario_base,
                                               delta0,
                                               nsim,
                                               alpha,
                                               seed,
                                               parallel,
                                               ncores,
                                               robust,
                                               eps,
                                               gamma_li,
                                               gate_c,
                                               gate_tau,
                                               gamma_mcp,
                                               delta_bounds,
                                               n_grid_opt,
                                               cl = NULL) {
  sc <- scenario_base
  sc$theta0 <- 0
  sc$delta0 <- delta0
  zcrit <- stats::qnorm(alpha)

  sim_one <- function(s) {
    sim <- simulate_hybrid_cox(
      nI1 = sc$nI1, nI0 = sc$nI0, nE = sc$nE,
      theta0 = sc$theta0, delta0 = sc$delta0,
      p = sc$p, beta = sc$beta, rho = sc$rho,
      cov_shift = sc$cov_shift, shape = sc$shape,
      lambda = sc$lambda, target_cens = sc$target_cens
    )
    .calibration_rejects_for_dataset(
      dat = sim$data,
      method = method,
      lambda_grid = lambda_grid,
      zcrit = zcrit,
      gamma_li = gamma_li,
      gate_c = gate_c,
      gate_tau = gate_tau,
      gamma_mcp = gamma_mcp,
      delta_bounds = delta_bounds,
      robust = robust,
      eps = eps,
      n_grid_opt = n_grid_opt
    )
  }

  if (!is.null(seed)) {
    if (!parallel) set.seed(seed)
  }

  if (parallel) {
    if (is.null(cl)) {
      stop("Internal error: parallel=TRUE requires a cluster object")
    }
    vals <- parallel::parLapply(cl, seq_len(nsim), sim_one)
  } else {
    vals <- lapply(seq_len(nsim), sim_one)
  }

  raw <- do.call(rbind, vals)
  if (is.null(raw) || nrow(raw) == 0) {
    return(data.frame(
      method = method, lambda = lambda_grid,
      delta0 = delta0, driftHR = exp(delta0),
      type1_model_based = NA_real_, type1_sandwich = NA_real_,
      n_nonmissing_model_based = 0, n_nonmissing_sandwich = 0
    ))
  }

  out <- lapply(lambda_grid, function(lam) {
    sub <- raw[abs(raw$lambda - lam) <= max(1e-12, abs(lam) * 1e-12),
               , drop = FALSE]
    data.frame(
      method = method,
      lambda = lam,
      delta0 = delta0,
      driftHR = exp(delta0),
      type1_model_based = mean(sub$reject_mod, na.rm = TRUE),
      type1_sandwich    = mean(sub$reject_sand, na.rm = TRUE),
      n_nonmissing_model_based = sum(is.finite(sub$reject_mod)),
      n_nonmissing_sandwich    = sum(is.finite(sub$reject_sand)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Build the worst-drift summary table across a lambda grid
#'
#' @keywords internal
#' @noRd
.summarize_calibration_details <- function(details_tbl, nsim) {
  if (is.null(details_tbl) || nrow(details_tbl) == 0) return(data.frame())

  lambda_vals <- sort(unique(details_tbl$lambda))
  out <- lapply(lambda_vals, function(lam) {
    sub <- details_tbl[abs(details_tbl$lambda - lam) <=
                         max(1e-12, abs(lam) * 1e-12), , drop = FALSE]
    worst_mod  <- max(sub$type1_model_based, na.rm = TRUE)
    worst_sand <- max(sub$type1_sandwich,    na.rm = TRUE)
    data.frame(
      method = unique(sub$method)[1],
      lambda = lam,
      worst_type1_model_based = worst_mod,
      worst_type1_sandwich    = worst_sand,
      worst_type1_model_based_upper95 = .calibration_mc_upper(worst_mod, nsim),
      worst_type1_sandwich_upper95    = .calibration_mc_upper(worst_sand, nsim),
      n_drift_evaluated = length(unique(sub$delta0)),
      nsim = nsim,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Convert a calibration summary table to long format
#'
#' @keywords internal
#' @noRd
.calibration_table_long <- function(method,
                                    summary_tbl,
                                    alpha,
                                    alpha_cal,
                                    lambda_star_model_based,
                                    lambda_star_sandwich) {
  if (is.null(summary_tbl) || nrow(summary_tbl) == 0) return(data.frame())
  rbind(
    data.frame(
      method = method, inference = "model_based",
      lambda = summary_tbl$lambda,
      worst_type1         = summary_tbl$worst_type1_model_based,
      worst_type1_upper95 = summary_tbl$worst_type1_model_based_upper95,
      lambda_star         = lambda_star_model_based,
      alpha = alpha, alpha_cal = alpha_cal,
      stage = summary_tbl$stage,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = method, inference = "sandwich",
      lambda = summary_tbl$lambda,
      worst_type1         = summary_tbl$worst_type1_sandwich,
      worst_type1_upper95 = summary_tbl$worst_type1_sandwich_upper95,
      lambda_star         = lambda_star_sandwich,
      alpha = alpha, alpha_cal = alpha_cal,
      stage = summary_tbl$stage,
      stringsAsFactors = FALSE
    )
  )
}

#' Start a parallel cluster for calibration
#'
#' @keywords internal
#' @noRd
.start_calibration_cluster <- function(ncores = NULL, seed = NULL) {
  if (is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl, library(fdb))
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
  cl
}

# User-facing calibration functions =====================================

#' Single-stage lambda calibration for one borrowing method
#'
#' For each lambda in \code{lambda_grid}, simulates \code{nsim}
#' replicates under each drift value in \code{drift_set} (with
#' \eqn{\theta_0 = 0}), records the model-based and sandwich rejection
#' rates, and selects the largest lambda whose worst-case rejection
#' rate over the drift set does not exceed \code{alpha_cal}.
#'
#' @param method One of \code{"Li"}, \code{"P1"}, \code{"P2"},
#'   \code{"P3"}, \code{"P4"}.
#' @param lambda_grid Numeric vector of candidate lambda values.
#' @param scenario_base A scenario list as accepted by
#'   \code{\link{run_simulation}}. Its \code{theta0} and \code{delta0}
#'   are overwritten internally.
#' @param drift_set Numeric vector of drift values (log HR) on which
#'   to evaluate type I error.
#' @param nsim Number of replicates per drift value.
#' @param alpha Nominal level used for rejection decisions.
#' @param alpha_cal Calibration threshold (worst-case rejection rate
#'   must not exceed this).
#' @param seed RNG seed.
#' @param parallel Logical; enable parallel evaluation across replicates.
#' @param ncores Number of cores for parallel evaluation.
#' @param robust Use robust (Lin-Wei) Cox SEs.
#' @param eps Smoothing parameter.
#' @param gamma_li Adaptive lasso exponent.
#' @param gate_c,gate_tau P2 gate parameters.
#' @param gamma_mcp MCP shape parameter for P3.
#' @param delta_bounds Optimization interval for delta.
#' @param n_grid_opt Coarse-grid points for non-convex objectives.
#' @param early_stop_drift Logical; if \code{TRUE}, drop lambda values
#'   that have already failed the calibration constraint after
#'   evaluating a subset of drift values.
#' @param stop_rule Either \code{"point"} (use point estimate of
#'   worst-case rejection rate for early stopping) or \code{"upper95"}
#'   (use the Monte Carlo 95\% upper confidence bound).
#' @param select_rule Selection rule for the final calibrated lambda
#'   (same options as \code{stop_rule}).
#' @return A list with elements \code{method}, \code{details}
#'   (per-lambda x drift), \code{summary} (per-lambda worst-case),
#'   \code{calibration_table} (long format), and \code{lambda_star}
#'   (selected lambdas for both inference types).
#'
#' @examples
#' \donttest{
#' cal <- calibrate_lambda_grid(method = "P1",
#'                              lambda_grid = exp(seq(log(0.05), log(1), length.out = 4)),
#'                              scenario_base = scenario_S1,
#'                              drift_set = make_drift_set_from_values(
#'                                c(0.9, 1.0, 1.1)),
#'                              nsim = 100, seed = 1)
#' cal$lambda_star
#' }
#' @export
calibrate_lambda_grid <- function(method = c("Li", "P1", "P2", "P3", "P4"),
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
                                  select_rule = c("point", "upper95")) {

  method <- match.arg(method)
  stop_rule <- match.arg(stop_rule)
  select_rule <- match.arg(select_rule)
  stopifnot(is.numeric(lambda_grid), length(lambda_grid) >= 1,
            nsim > 0, alpha > 0, alpha < 1,
            alpha_cal > 0, alpha_cal < 1)

  lambda_grid <- sort(unique(lambda_grid))
  drift_set <- sort(unique(drift_set))

  cl <- NULL
  if (parallel) {
    cl <- .start_calibration_cluster(ncores = ncores, seed = seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else if (!is.null(seed)) {
    set.seed(seed)
  }

  active_lambdas <- lambda_grid
  drift_details <- list()
  detail_id <- 0

  for (j in seq_along(drift_set)) {
    if (length(active_lambdas) == 0) break

    drift_seed <- if (is.null(seed)) NULL else seed + 10000L * j
    drift_tbl <- .evaluate_calibration_drift_cached(
      method = method,
      lambda_grid = active_lambdas,
      scenario_base = scenario_base,
      delta0 = drift_set[j],
      nsim = nsim,
      alpha = alpha,
      seed = drift_seed,
      parallel = parallel,
      ncores = ncores,
      robust = robust,
      eps = eps,
      gamma_li = gamma_li,
      gate_c = gate_c,
      gate_tau = gate_tau,
      gamma_mcp = gamma_mcp,
      delta_bounds = delta_bounds,
      n_grid_opt = n_grid_opt,
      cl = cl
    )

    detail_id <- detail_id + 1
    drift_details[[detail_id]] <- drift_tbl

    if (early_stop_drift) {
      current_details <- do.call(rbind, drift_details)
      current_summary <- .summarize_calibration_details(current_details,
                                                        nsim = nsim)

      if (stop_rule == "upper95") {
        keep_mod  <- current_summary$worst_type1_model_based_upper95 <= alpha_cal
        keep_sand <- current_summary$worst_type1_sandwich_upper95    <= alpha_cal
      } else {
        keep_mod  <- current_summary$worst_type1_model_based <= alpha_cal
        keep_sand <- current_summary$worst_type1_sandwich    <= alpha_cal
      }

      keep <- keep_mod | keep_sand
      active_lambdas <- current_summary$lambda[
        keep & current_summary$lambda %in% active_lambdas
      ]
      if (length(active_lambdas) == 0) {
        message(sprintf(
          "Drift-level early stopping: all lambdas failed for method=%s",
          method))
      }
    }
  }

  details_tbl <- if (length(drift_details) == 0) {
    data.frame()
  } else {
    do.call(rbind, drift_details)
  }
  summary_tbl <- .summarize_calibration_details(details_tbl, nsim = nsim)

  lambda_star_model_based <- .select_lambda_star(summary_tbl, alpha_cal,
                                                 "model_based", select_rule)
  lambda_star_sandwich    <- .select_lambda_star(summary_tbl, alpha_cal,
                                                 "sandwich",    select_rule)

  summary_tbl$lambda_star_model_based <- lambda_star_model_based
  summary_tbl$lambda_star_sandwich    <- lambda_star_sandwich
  summary_tbl$alpha <- alpha
  summary_tbl$alpha_cal <- alpha_cal
  summary_tbl$selection_rule <- select_rule
  summary_tbl$stage <- "single"

  calibration_table_long <- .calibration_table_long(
    method = method,
    summary_tbl = summary_tbl,
    alpha = alpha, alpha_cal = alpha_cal,
    lambda_star_model_based = lambda_star_model_based,
    lambda_star_sandwich    = lambda_star_sandwich
  )

  list(
    method  = method,
    details = details_tbl,
    summary = summary_tbl,
    calibration_table = calibration_table_long,
    lambda_star = data.frame(
      method = method,
      inference = c("model_based", "sandwich"),
      lambda_star = c(lambda_star_model_based, lambda_star_sandwich),
      stringsAsFactors = FALSE
    )
  )
}

#' Two-stage lambda calibration for one borrowing method
#'
#' Performs a coarse-grid calibration on a reduced drift set
#' (\code{drift_set_cal}), refines the grid around the calibrated
#' lambda, repeats the fine-grid calibration on \code{drift_set_cal},
#' and optionally confirms the final lambda on a wider drift set
#' \code{drift_set_confirm}.
#'
#' @param method One of \code{"Li"}, \code{"P1"}, \code{"P2"},
#'   \code{"P3"}, \code{"P4"}.
#' @param lambda_grid_coarse Initial coarse lambda grid.
#' @param scenario_base A scenario list.
#' @param drift_set_cal Reduced drift set for calibration.
#' @param drift_set_confirm Wider drift set for confirmation
#'   (defaults to \code{NULL}, which disables confirmation).
#' @param nsim_cal Replicates per drift value in the calibration stages.
#' @param nsim_confirm Replicates per drift value in the confirmation
#'   stage.
#' @param alpha,alpha_cal Nominal level and calibration threshold.
#' @param seed RNG seed (subsequent stages use offsets of this seed).
#' @param parallel,ncores Parallelization controls.
#' @param robust,eps,gamma_li,gate_c,gate_tau,gamma_mcp,delta_bounds,n_grid_opt
#'   Same as in \code{\link{calibrate_lambda_grid}}.
#' @param n_fine Number of points in the fine lambda grid.
#' @param primary_inference Which inference type drives the
#'   refinement (default \code{"sandwich"}).
#' @param confirm_full_drift Logical; perform the confirmation stage
#'   on \code{drift_set_confirm}.
#' @param early_stop_drift,stop_rule,select_rule As in
#'   \code{\link{calibrate_lambda_grid}}.
#' @return A list with components \code{method}, \code{coarse},
#'   \code{fine}, \code{confirm}, \code{calibration_summary},
#'   \code{calibration_table}, \code{lambda_star}, and
#'   \code{final_stage}.
#' @export
calibrate_lambda_grid_two_stage <- function(method = c("Li", "P1", "P2", "P3", "P4"),
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
                                            select_rule = c("point", "upper95")) {
  method <- match.arg(method)
  primary_inference <- match.arg(primary_inference)
  stop_rule <- match.arg(stop_rule)
  select_rule <- match.arg(select_rule)

  stage1 <- calibrate_lambda_grid(
    method = method, lambda_grid = lambda_grid_coarse,
    scenario_base = scenario_base, drift_set = drift_set_cal,
    nsim = nsim_cal, alpha = alpha, alpha_cal = alpha_cal,
    seed = seed + 11, parallel = parallel, ncores = ncores,
    robust = robust, eps = eps,
    gamma_li = gamma_li, gate_c = gate_c, gate_tau = gate_tau,
    gamma_mcp = gamma_mcp, delta_bounds = delta_bounds,
    n_grid_opt = n_grid_opt, early_stop_drift = early_stop_drift,
    stop_rule = stop_rule, select_rule = select_rule
  )
  stage1$summary$stage <- "coarse_reduced_drift"
  stage1$calibration_table$stage <- "coarse_reduced_drift"

  primary_star1 <- stage1$lambda_star$lambda_star[
    stage1$lambda_star$inference == primary_inference
  ]
  lambda_grid_fine <- .make_refined_lambda_grid(lambda_grid_coarse,
                                                primary_star1,
                                                n_fine = n_fine)

  stage2 <- calibrate_lambda_grid(
    method = method, lambda_grid = lambda_grid_fine,
    scenario_base = scenario_base, drift_set = drift_set_cal,
    nsim = nsim_cal, alpha = alpha, alpha_cal = alpha_cal,
    seed = seed + 22, parallel = parallel, ncores = ncores,
    robust = robust, eps = eps,
    gamma_li = gamma_li, gate_c = gate_c, gate_tau = gate_tau,
    gamma_mcp = gamma_mcp, delta_bounds = delta_bounds,
    n_grid_opt = n_grid_opt, early_stop_drift = early_stop_drift,
    stop_rule = stop_rule, select_rule = select_rule
  )
  stage2$summary$stage <- "fine_reduced_drift"
  stage2$calibration_table$stage <- "fine_reduced_drift"

  final_stage <- stage2
  confirm <- NULL
  if (isTRUE(confirm_full_drift) && !is.null(drift_set_confirm)) {
    primary_star2 <- stage2$lambda_star$lambda_star[
      stage2$lambda_star$inference == primary_inference
    ]
    if (is.finite(primary_star2)) {
      lambda_confirm <- sort(unique(
        lambda_grid_fine[lambda_grid_fine <= primary_star2 * (1 + 1e-12)]
      ))
      if (length(lambda_confirm) == 0) lambda_confirm <- primary_star2
    } else {
      lambda_confirm <- sort(unique(lambda_grid_fine))
    }

    confirm <- calibrate_lambda_grid(
      method = method, lambda_grid = lambda_confirm,
      scenario_base = scenario_base, drift_set = drift_set_confirm,
      nsim = nsim_confirm, alpha = alpha, alpha_cal = alpha_cal,
      seed = seed + 33, parallel = parallel, ncores = ncores,
      robust = robust, eps = eps,
      gamma_li = gamma_li, gate_c = gate_c, gate_tau = gate_tau,
      gamma_mcp = gamma_mcp, delta_bounds = delta_bounds,
      n_grid_opt = n_grid_opt, early_stop_drift = early_stop_drift,
      stop_rule = stop_rule, select_rule = select_rule
    )
    confirm$summary$stage <- "confirm_full_drift"
    confirm$calibration_table$stage <- "confirm_full_drift"
    final_stage <- confirm
  }

  combined_summary <- do.call(
    rbind,
    Filter(Negate(is.null),
           list(stage1$summary, stage2$summary,
                if (!is.null(confirm)) confirm$summary else NULL))
  )
  combined_table <- do.call(
    rbind,
    Filter(Negate(is.null),
           list(stage1$calibration_table, stage2$calibration_table,
                if (!is.null(confirm)) confirm$calibration_table else NULL))
  )

  list(method = method,
       coarse = stage1, fine = stage2, confirm = confirm,
       calibration_summary = combined_summary,
       calibration_table   = combined_table,
       lambda_star         = final_stage$lambda_star,
       final_stage         = final_stage)
}

#' Calibrate lambda for all borrowing methods
#'
#' Runs lambda calibration for all five penalized borrowing methods
#' (Li adaptive lasso plus P1-P4) and returns a unified set of
#' calibrated lambdas under both model-based and sandwich inference.
#'
#' @param lambda_grid Initial (coarse) lambda grid used for all
#'   methods.
#' @param scenario_base A scenario list.
#' @param drift_set The full drift set on which to evaluate the final
#'   confirmation (also used as the default calibration drift set).
#' @param drift_set_cal Reduced drift set for the calibration stages.
#' @param drift_set_confirm Drift set for the optional confirmation
#'   stage (defaults to \code{drift_set}).
#' @param nsim_cal,nsim_confirm Replicates per drift value in
#'   calibration and confirmation stages.
#' @param alpha,alpha_cal,seed,parallel,ncores,robust,eps,gamma_li,
#'   gate_c,gate_tau,gamma_mcp,delta_bounds,n_grid_opt,n_fine,
#'   primary_inference,confirm_full_drift,early_stop_drift,stop_rule,
#'   select_rule As in \code{\link{calibrate_lambda_grid_two_stage}}.
#' @param two_stage Logical; use two-stage calibration if \code{TRUE},
#'   otherwise single-stage.
#' @return A list with \code{calibration_table} (long form),
#'   \code{calibration_summary} (per-method, per-lambda worst case),
#'   \code{lambda_star} (5 x 2 matrix indexed by method and
#'   inference), \code{lambda_star_table} (long form), and
#'   \code{method_calibrations} (per-method full output).
#' @export
calibrate_all_lambdas <- function(lambda_grid,
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
                                  select_rule = c("point", "upper95")) {
  primary_inference <- match.arg(primary_inference)
  stop_rule <- match.arg(stop_rule)
  select_rule <- match.arg(select_rule)

  if (is.null(drift_set_cal)) drift_set_cal <- drift_set
  methods <- c("Li", "P1", "P2", "P3", "P4")

  cal_list <- vector("list", length(methods))
  for (k in seq_along(methods)) {
    message(sprintf("Starting calibration for method=%s", methods[k]))
    if (isTRUE(two_stage)) {
      cal_list[[k]] <- calibrate_lambda_grid_two_stage(
        method = methods[k],
        lambda_grid_coarse = lambda_grid,
        scenario_base = scenario_base,
        drift_set_cal = drift_set_cal,
        drift_set_confirm = drift_set_confirm,
        nsim_cal = nsim_cal, nsim_confirm = nsim_confirm,
        alpha = alpha, alpha_cal = alpha_cal,
        seed = seed + k * 1000,
        parallel = parallel, ncores = ncores,
        robust = robust, eps = eps,
        gamma_li = gamma_li, gate_c = gate_c, gate_tau = gate_tau,
        gamma_mcp = gamma_mcp, delta_bounds = delta_bounds,
        n_grid_opt = n_grid_opt, n_fine = n_fine,
        primary_inference = primary_inference,
        confirm_full_drift = confirm_full_drift,
        early_stop_drift = early_stop_drift,
        stop_rule = stop_rule, select_rule = select_rule
      )
    } else {
      one <- calibrate_lambda_grid(
        method = methods[k], lambda_grid = lambda_grid,
        scenario_base = scenario_base, drift_set = drift_set_cal,
        nsim = nsim_cal, alpha = alpha, alpha_cal = alpha_cal,
        seed = seed + k * 1000,
        parallel = parallel, ncores = ncores,
        robust = robust, eps = eps,
        gamma_li = gamma_li, gate_c = gate_c, gate_tau = gate_tau,
        gamma_mcp = gamma_mcp, delta_bounds = delta_bounds,
        n_grid_opt = n_grid_opt,
        early_stop_drift = early_stop_drift,
        stop_rule = stop_rule, select_rule = select_rule
      )
      cal_list[[k]] <- list(
        method = methods[k],
        calibration_summary = one$summary,
        calibration_table = one$calibration_table,
        lambda_star = one$lambda_star, final_stage = one
      )
    }
  }

  cal_tbl     <- do.call(rbind, lapply(cal_list, function(x) x$calibration_table))
  cal_summary <- do.call(rbind, lapply(cal_list, function(x) x$calibration_summary))
  star_tbl    <- do.call(rbind, lapply(cal_list, function(x) x$lambda_star))

  lambda_star <- matrix(NA_real_, nrow = length(methods), ncol = 2,
                        dimnames = list(methods,
                                        c("model_based", "sandwich")))
  for (i in seq_len(nrow(star_tbl))) {
    lambda_star[star_tbl$method[i], star_tbl$inference[i]] <-
      star_tbl$lambda_star[i]
  }

  list(
    calibration_table   = cal_tbl,
    calibration_summary = cal_summary,
    lambda_star         = lambda_star,
    lambda_star_table   = star_tbl,
    method_calibrations = cal_list
  )
}
