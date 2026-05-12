# Profile-likelihood optimization for delta ==============================

#' One-dimensional minimization of a penalized profile objective in delta
#'
#' For convex or near-convex objectives, \code{n_grid <= 1} uses
#' \code{stats::optimize} over the full interval. For nonconvex or
#' evidence-gated penalties (notably the gated L1 (P2) and MCP (P3)
#' penalties), a coarse grid is evaluated first, the best grid point
#' is chosen, and \code{stats::optimize} is then run over the
#' neighboring interval.
#'
#' @param dat A data frame with columns \code{time}, \code{status},
#'   \code{T}, \code{Z}, and covariates named in \code{xnames}.
#' @param xnames Character vector of covariate names.
#' @param obj_fn Objective function of a single argument \code{delta}.
#' @param delta_bounds Numeric vector of length 2, \code{(lower, upper)}.
#' @param robust If \code{TRUE}, use robust SEs from \code{coxph}.
#' @param n_grid Number of grid points for the coarse search. Use
#'   \code{n_grid <= 1} to disable grid-assisted search.
#' @return A list with the optimized \code{delta_hat}, profiled
#'   estimates \code{theta_hat}, \code{se_theta}, \code{beta_hat},
#'   the profile log-likelihood \code{loglik}, the optimum objective
#'   value \code{opt_obj}, and the number of grid points used.
#' @keywords internal
#' @noRd
optimize_delta <- function(dat, xnames, obj_fn,
                           delta_bounds = DEFAULT_DELTA_BOUNDS,
                           robust = FALSE,
                           n_grid = 1) {
  if (length(delta_bounds) != 2 || !all(is.finite(delta_bounds)) ||
      delta_bounds[1] >= delta_bounds[2]) {
    stop("delta_bounds must be a finite increasing vector of length 2")
  }

  n_grid <- as.integer(n_grid)

  safe_obj <- function(d) {
    val <- tryCatch(obj_fn(d), error = function(e) Inf)
    if (!is.finite(val)) Inf else val
  }

  if (is.na(n_grid) || n_grid <= 1) {
    opt <- stats::optimize(safe_obj, interval = delta_bounds)
  } else {
    n_grid <- max(3L, n_grid)
    grid <- seq(delta_bounds[1], delta_bounds[2], length.out = n_grid)
    vals <- vapply(grid, safe_obj, numeric(1))

    finite_idx <- which(is.finite(vals))
    if (length(finite_idx) == 0) {
      stop("All objective evaluations failed in grid-assisted optimize_delta().")
    }

    best_idx <- finite_idx[which.min(vals[finite_idx])]
    lo_idx <- max(1L, best_idx - 1L)
    hi_idx <- min(length(grid), best_idx + 1L)
    local_interval <- c(grid[lo_idx], grid[hi_idx])

    if (local_interval[1] == local_interval[2]) {
      opt <- list(minimum = grid[best_idx], objective = vals[best_idx])
    } else {
      opt <- stats::optimize(safe_obj, interval = local_interval)
      if (!is.finite(opt$objective) || vals[best_idx] < opt$objective) {
        opt <- list(minimum = grid[best_idx], objective = vals[best_idx])
      }
    }
  }

  delta_hat <- opt$minimum
  prof <- cox_fit_profile(dat, delta_hat, xnames, robust = robust)
  list(
    delta_hat  = delta_hat,
    theta_hat  = prof$theta_hat,
    se_theta   = prof$se_theta,
    beta_hat   = prof$beta_hat,
    loglik     = prof$loglik,
    opt_obj    = opt$objective,
    n_grid_opt = n_grid
  )
}

#' Assemble a method result row including sandwich SE
#'
#' Wraps a method's profiled solution with sandwich variance estimation
#' and returns a list of estimates and inference quantities. The
#' \code{pen_curv} input is the penalty curvature
#' \eqn{p''_{\lambda,\varepsilon}(\hat\delta)} entering the sandwich bread.
#'
#' @keywords internal
#' @noRd
build_method_result <- function(method_name, sol, dat, xnames, pen_curv,
                                extra = list()) {
  z_mod <- sol$theta_hat / sol$se_theta

  sw <- tryCatch(
    compute_sandwich_se(dat, sol$delta_hat, sol$theta_hat,
                        sol$beta_hat, pen_curv, xnames),
    error = function(e) list(se_sand = NA_real_, note = conditionMessage(e))
  )
  se_sand <- sw$se_sand
  z_sand  <- if (is.finite(se_sand) && se_sand > 0) {
    sol$theta_hat / se_sand
  } else {
    NA_real_
  }

  c(
    list(method    = method_name,
         theta_hat = sol$theta_hat,
         se_theta  = sol$se_theta,
         se_sand   = se_sand,
         z         = z_mod,
         z_sand    = z_sand,
         delta_hat = sol$delta_hat,
         pen_curv  = pen_curv),
    extra
  )
}
