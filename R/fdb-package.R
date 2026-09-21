#' fdb: Frequentist Dynamic Borrowing for Hybrid-Control Survival Trials
#'
#' Implements likelihood-informed frequentist dynamic borrowing methods
#' for hybrid-control survival trials based on penalized Cox partial
#' likelihood estimation, together with design-stage lambda calibration
#' and a simulation harness for evaluating operating characteristics.
#'
#' @section Main user-facing functions:
#' \itemize{
#'   \item \code{\link{simulate_hybrid_cox}}: Generate a hybrid-control
#'     Cox proportional hazards dataset.
#'   \item \code{\link{fit_all_methods}}: Fit all borrowing methods on a
#'     single dataset, returning both model-based and sandwich-based
#'     inference quantities.
#'   \item \code{\link{fit_one_penalized_method}}: Fit a single
#'     penalized borrowing method.
#'   \item \code{\link{run_simulation}}: Run a Monte Carlo simulation
#'     under a fixed scenario, comparing all methods.
#'   \item \code{\link{calibrate_lambda_grid}},
#'     \code{\link{calibrate_lambda_grid_two_stage}},
#'     \code{\link{calibrate_all_lambdas}}: Design-stage lambda
#'     calibration utilities.
#'   \item \code{\link{run_fdb_study}}: One-stop wrapper that performs
#'     lambda calibration and then evaluates type I error and power
#'     curves across population drift scenarios.
#' }
#'
#' @section Penalty methods:
#' The package implements four likelihood-informed penalties together
#' with the adaptive lasso comparator of Li et al. (2023):
#' \itemize{
#'   \item \strong{LiAdaptiveLasso}: adaptive lasso (Li et al. 2023).
#'   \item \strong{P1_SEScaledL1}: precision-weighted L1 penalty.
#'   \item \strong{P2_GatedL1}: smoothed integrated-gate penalty.
#'   \item \strong{P3_SEScaledMCP}: information-adaptive minimax
#'     concave penalty (MCP).
#'   \item \strong{P4_LRWeightedL1}: likelihood-ratio-weighted L1
#'     penalty.
#' }
#'
#' @section References:
#' Li, R., Lin, R., Huang, J., Tian, L., and Zhu, J. (2023).
#' A frequentist approach to dynamic borrowing.
#' \emph{Biometrical Journal} 65(7), 2100406.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under
#' minimax concave penalty. \emph{The Annals of Statistics} 38(2),
#' 894-942.
#'
#' Andersen, P. K. and Gill, R. D. (1982). Cox's regression model for
#' counting processes: A large sample study. \emph{The Annals of
#' Statistics} 10(4), 1100-1120.
#'
#' @docType package
#' @name fdb-package
#' @aliases fdb
NULL

# Global numerical constants ---------------------------------------------

#' @keywords internal
#' @noRd
NUMERIC_ZERO <- 1e-10

#' @keywords internal
#' @noRd
DEFAULT_DELTA_BOUNDS <- c(-3, 3)

#' @keywords internal
#' @noRd
SMOOTH_EPS <- 1e-3

# Software default; choose and report this value before calibration.
DEFAULT_RHO_MCP <- 0.1

#' @keywords internal
#' @noRd
DEFAULT_N_GRID_OPT <- 21

# Null-coalescing utility -------------------------------------------------

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
