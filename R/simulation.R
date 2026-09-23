# Simulation driver and ESS computation ==================================

#' Run a Monte Carlo simulation under a fixed scenario
#'
#' For each of \code{nsim} replicates, simulates a hybrid-control Cox
#' dataset and fits all borrowing methods. Returns per-replicate raw
#' results and aggregated summaries under both model-based and
#' sandwich inference.
#'
#' @param nsim Number of Monte Carlo replicates, at least two.
#' @param scenario A list of data-generating parameters with elements
#'   \code{nI1}, \code{nI0}, \code{nE}, \code{theta0}, \code{delta0},
#'   \code{p}, \code{beta}, \code{rho}, \code{cov_shift},
#'   \code{shape}, \code{lambda}, and \code{target_cens}.
#'   See \code{\link{scenario_S1}} for an example.
#' @param lambdas A list of tuning parameters with elements
#'   \code{lambda_li}, \code{gamma_li}, \code{lambda_p1},
#'   \code{lambda_p2}, \code{gate_c}, \code{gate_tau},
#'   \code{lambda_p3}, \code{gamma_mcp}, \code{rho_mcp} (optional, default 0.1), \code{lambda_p4}, and
#'   \code{delta_bounds}. See \code{\link{lambdas_default}} for an
#'   example.
#' @param alpha Nominal one-sided (or two-sided) significance level.
#' @param one_sided Logical; if \code{TRUE}, rejection corresponds to
#'   a hazard ratio less than one.
#' @param seed RNG seed.
#' @param parallel Logical; if \code{TRUE}, parallelize replicates
#'   using \code{parallel::parLapply}.
#' @param ncores Number of workers; \code{NULL} uses two. Checks use at most two.
#' @param robust Use robust (Lin-Wei) Cox standard errors.
#' @param eps Smoothing parameter for \eqn{|\delta|_\varepsilon}.
#' @param n_grid_opt Number of grid points for the coarse search in
#'   non-convex objectives.
#' @return A list with:
#'   \describe{
#'     \item{\code{raw}}{Per-replicate results data frame.}
#'     \item{\code{summary}}{Aggregated summaries for both inference
#'       types (valid-result counts, Monte Carlo standard errors, rejection rate,
#'       bias, RMSE, empirical SE, average SE,
#'       95\% coverage), with an \code{inference} column.}
#'     \item{\code{scenario}, \code{lambdas}, \code{settings}}{The
#'       inputs and run metadata.}
#'   }
#'
#' @examples
#' \donttest{
#' sim_out <- run_simulation(nsim = 2, scenario = scenario_S1,
#'                           lambdas = lambdas_default, alpha = 0.025, seed = 1)
#' subset(sim_out$summary, inference == "sandwich")
#' }
#' @export
run_simulation <- function(nsim = 200,
                           scenario,
                           lambdas,
                           alpha = 0.025,
                           one_sided = TRUE,
                           seed = 1,
                           parallel = FALSE,
                           ncores = NULL,
                           robust = FALSE,
                           eps = SMOOTH_EPS,
                           n_grid_opt = DEFAULT_N_GRID_OPT) {

  stopifnot(is.list(scenario), is.list(lambdas),
            nsim > 0, alpha > 0, alpha < 1)

  .validate_count(nsim, "nsim", 2L)

  if (!is.null(seed)) {
    if (parallel) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    } else {
      set.seed(seed)
    }
  }

  # Evaluate caller settings before sending the closure to PSOCK workers.
  force(robust); force(eps); force(n_grid_opt); force(one_sided)
  zcrit <- stats::qnorm(alpha)

  sim_one_rep <- function(s) {
    sim <- simulate_hybrid_cox(
      nI1 = scenario$nI1, nI0 = scenario$nI0, nE = scenario$nE,
      theta0 = scenario$theta0, delta0 = scenario$delta0,
      p = scenario$p, beta = scenario$beta, rho = scenario$rho,
      cov_shift = scenario$cov_shift, shape = scenario$shape,
      lambda = scenario$lambda, target_cens = scenario$target_cens
    )
    dat <- sim$data

    res <- fit_all_methods(
      dat,
      lambda_li = lambdas$lambda_li, gamma_li = lambdas$gamma_li,
      lambda_p1 = lambdas$lambda_p1,
      lambda_p2 = lambdas$lambda_p2, gate_c = lambdas$gate_c,
      gate_tau = lambdas$gate_tau,
      lambda_p3 = lambdas$lambda_p3, gamma_mcp = lambdas$gamma_mcp,
      rho_mcp = lambdas$rho_mcp %||% DEFAULT_RHO_MCP,
      lambda_p4 = lambdas$lambda_p4,
      delta_bounds = lambdas$delta_bounds %||% DEFAULT_DELTA_BOUNDS,
      robust = robust, eps = eps, n_grid_opt = n_grid_opt
    )

    res$reject_mod  <- if (one_sided) (res$z < zcrit) else (abs(res$z) > stats::qnorm(1 - alpha / 2))
    res$reject_sand <- if (one_sided) (res$z_sand < zcrit) else (abs(res$z_sand) > stats::qnorm(1 - alpha / 2))

    res$sim    <- s
    res$theta0 <- scenario$theta0
    res$delta0 <- scenario$delta0
    res
  }

  if (parallel) {
    ncores <- .resolve_ncores(ncores)
    cl <- .start_calibration_cluster(ncores, seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    all_res <- parallel::parLapply(cl, seq_len(nsim), sim_one_rep)
  } else {
    all_res <- lapply(seq_len(nsim), sim_one_rep)
  }

  out <- do.call(rbind, all_res)

  summ <- rbind(.simulation_summary(out, scenario$theta0, "model_based"),
                .simulation_summary(out, scenario$theta0, "sandwich"))
  if (any(summ$n_missing > 0)) {
    warning("Some simulation results are missing; inspect n_valid and n_missing.", call. = FALSE)
  }

  list(
    raw      = out,
    summary  = summ,
    scenario = scenario,
    lambdas  = lambdas,
    settings = list(nsim = nsim, alpha = alpha, one_sided = one_sided,
                    seed = seed, parallel = parallel, robust = robust,
                    eps = eps, n_grid_opt = n_grid_opt,
                    timestamp = Sys.time(),
                    session_info = utils::sessionInfo())
  )
}

# Effective sample size ==================================================

#' Compute effective sample size from raw simulation output
#'
#' Computes a variance-ratio effective sample size for each method
#' relative to a reference method (typically internal-only), based on
#' the Monte Carlo variance of the treatment effect estimator across
#' replicates: \eqn{ESS_m = NS \cdot (Var_{ref} / Var_m - 1)}.
#'
#' @param raw_df A data frame with at least \code{method} and
#'   \code{theta_hat} columns.
#' @param NS Reference sample-size scale. The study wrappers use the total
#'   internal randomized sample size, \code{nI1 + nI0}.
#' @param ref_method Name of the reference method (default
#'   \code{"InternalOnly"}).
#' @param methods_exclude Optional character vector of methods to
#'   exclude from the output.
#' @return A data frame with \code{method} and \code{ESS} columns. ESS is a
#'   variance-equivalent gain, can be negative, and does not measure bias.
#' @examples
#' \donttest{
#' sim_out <- run_simulation(nsim = 2, scenario = scenario_S1,
#'                           lambdas = lambdas_default, alpha = 0.025, seed = 1)
#' compute_ess_from_raw(sim_out$raw, NS = 300)
#' }
#' @export
compute_ess_from_raw <- function(raw_df,
                                 NS,
                                 ref_method = "InternalOnly",
                                 methods_exclude = NULL) {
  stopifnot(all(c("method", "theta_hat") %in% names(raw_df)))
  if (!ref_method %in% raw_df$method) {
    stop("Reference method not found: ", ref_method)
  }

  var_ref <- stats::var(raw_df$theta_hat[raw_df$method == ref_method],
                        na.rm = TRUE)
  if (!is.finite(var_ref) || var_ref <= 0) {
    stop("Reference variance is non-positive or not finite.")
  }

  methods <- unique(raw_df$method)
  if (!is.null(methods_exclude)) methods <- setdiff(methods, methods_exclude)

  ess <- lapply(methods, function(m) {
    v <- stats::var(raw_df$theta_hat[raw_df$method == m], na.rm = TRUE)
    ess_m <- if (!is.finite(v) || v <= 0) {
      NA_real_
    } else if (m == ref_method) {
      0
    } else {
      NS * (var_ref / v - 1)
    }
    data.frame(method = m, ESS = ess_m, stringsAsFactors = FALSE)
  })
  do.call(rbind, ess)
}

#' Append ESS to a simulation summary
#'
#' @param simres Output of \code{\link{run_simulation}}.
#' @param NS Reference sample-size scale. The study wrappers use the total
#'   internal randomized sample size, \code{nI1 + nI0}.
#' @param ref_method Reference method (default \code{"InternalOnly"}).
#' @return The same \code{simres} list with an \code{ESS} column
#'   merged into its \code{summary} component.
#' @export
add_ess_to_simulation_result <- function(simres, NS,
                                         ref_method = "InternalOnly") {
  ess_tbl <- compute_ess_from_raw(simres$raw, NS = NS,
                                  ref_method = ref_method)
  simres$summary <- merge(simres$summary, ess_tbl, by = "method", all.x = TRUE)
  simres
}
