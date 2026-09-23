# Internal controls shared by simulation and calibration.
.validate_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x)) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  invisible(x)
}

.offset_seed <- function(seed, offset) {
  if (is.null(seed)) NULL else seed + offset
}

.resolve_ncores <- function(ncores) {
  if (is.null(ncores)) ncores <- 2L
  .validate_count(ncores, "ncores")
  limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false"))
  if (limit %in% c("true", "warn") && ncores > 2L) {
    stop("Package checks must use at most two workers.", call. = FALSE)
  }
  as.integer(ncores)
}

.simulation_summary <- function(out, theta0, inference) {
  se_col <- if (inference == "sandwich") "se_sand" else "se_theta"
  reject_col <- if (inference == "sandwich") "reject_sand" else "reject_mod"
  mean_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  rows <- lapply(unique(out$method), function(method) {
    d <- out[out$method == method, , drop = FALSE]
    est_ok <- is.finite(d$theta_hat)
    se_ok <- is.finite(d[[se_col]]) & d[[se_col]] > 0
    valid <- est_ok & se_ok & is.finite(d[[reject_col]])
    errors <- d$theta_hat[est_ok] - theta0
    reject <- d[[reject_col]][valid]
    coverage <- abs(d$theta_hat[valid] - theta0) <= stats::qnorm(0.975) * d[[se_col]][valid]
    n_valid <- sum(valid)
    rate <- mean_finite(reject)
    cov <- mean_finite(as.numeric(coverage))
    data.frame(method = method, rej_rate = rate,
      bias = mean_finite(errors), rmse = sqrt(mean_finite(errors^2)),
      mse = mean_finite(errors^2), emp_se = stats::sd(d$theta_hat[est_ok]),
      avg_se = mean_finite(d[[se_col]][se_ok]), coverage_95 = cov,
      n_replicates = nrow(d), n_estimates = sum(est_ok), n_valid = n_valid,
      n_missing = nrow(d) - n_valid,
      mcse_rej_rate = if (n_valid) sqrt(rate * (1 - rate) / n_valid) else NA_real_,
      mcse_coverage = if (n_valid) sqrt(cov * (1 - cov) / n_valid) else NA_real_,
      inference = inference, stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, rows)
  ans[order(ans$method), , drop = FALSE]
}
