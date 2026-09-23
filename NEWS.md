# fdb 0.2.0

* Prepare source for CRAN checks and regenerate help from roxygen comments.
* Add optional `keep_raw` retention to curve and study workflows.
* Include valid/missing counts and Monte Carlo standard errors in simulation
  summaries; keep rows for methods with entirely missing inference results.
* Stop the study wrapper when calibration fails instead of substituting 0.20.
* Honor explicit `drift_set_confirm` grids. NULL defaults to the calibration
  grid. Existing calls with identical calibration/confirmation grids retain
  their calibration target; calls previously passing a different ignored
  confirmation grid now use that grid and can produce different tuning.
* Retain coarse candidates in the fine calibration search.
* Default parallel runs to two workers, propagate installed library paths,
  clean up failed cluster initialization, and support NULL seeds in wrappers.
* Document conditional model-based SEs, approximate sandwich inference,
  finite-grid calibration, residual undercoverage and signed ESS gains.
* Penalty formulas and fixed-lambda fitting algorithms are unchanged by this
  release-preparation pass. Earlier manuscript changes to P2/P3 still require
  recomputation of results from the original penalty definitions.

# fdb development version

* Replace P2 drift-times-gate with the manuscript integrated-gate penalty.
* Smooth the MCP origin and flat-tail transition; expose `rho_mcp` (default
  0.1) through fitting, calibration and study workflows.
* Match raw analytic curvature to the fitted penalties; retain clipping
  only in the plug-in variance calculation.
* Support no-covariate Cox fits and simulation (`p = 0`).
* Old P2/P3 fits and calibration results require recomputation.

# fdb 0.1.0

* Initial release.
* User-facing functions: `simulate_hybrid_cox`, `fit_internal_only`,
  `fit_naive_pooled`, `fit_li_adaptive_lasso`, `fit_P1_precision_L1`,
  `fit_P2_gated_L1`, `fit_P3_info_MCP`, `fit_P4_LRweighted_L1`,
  `fit_all_methods`, `fit_one_penalized_method`, `compute_sandwich_se`,
  `run_simulation`, `run_drift_curve`, `compute_ess_from_raw`,
  `add_ess_to_simulation_result`, `make_drift_set`,
  `make_drift_set_from_values`, `calibrate_lambda_grid`,
  `calibrate_lambda_grid_two_stage`, `calibrate_all_lambdas`,
  and `run_fdb_study`.
* Reference scenario `scenario_S1` and default tuning list
  `lambdas_default` exported for examples and tests.
