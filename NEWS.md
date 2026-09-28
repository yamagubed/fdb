# fdb 0.2.2


P3 now applies MCP to `delta / se_delta`, with lambda and gamma defined on
that standardized scale. The smooth penalty integrates its slope to
`(sqrt(delta^2 + eps^2) - eps) / se_delta`; `eps` remains in log-HR units.
The transition half-width is `rho_mcp * gamma_mcp * lambda` in SE units.
Both derivative chain factors are included in the raw curvature; negative
curvature is still clipped only for the sandwich bread. The optimizer includes
transition points and refines all detected local minima for P3.

This is a substantive replacement: recalibrate and reevaluate P3. Previous
P3 lambdas, estimates, and result files do not describe this method.
Other penalties and the fixed independent censoring design are unchanged.
The returned `lambda_eff` is retained as the origin slope `lambda / se_delta`;
it is no longer the lambda argument supplied to the MCP primitive.

# fdb 0.2.1

* Censoring now uses a fixed exponential rate calculated from the internal
  population distribution, independent of realized replicate event times.
  Added calibrate_censor_rate() and optional censor_rate. Target zero yields
  exactly no censoring. Removed the stochastic per-replicate rate search.
* Simulation and calibration pass the resolved fixed rate to every replicate.
  Existing simulation results and calibrated lambdas must not be mixed with
  results from this revised generator.
* Fit and simulation method labels now match the manuscript. Short calibration
  API keys Li/P1/P2/P3/P4 remain supported; method_labels() translates keys and
  legacy output labels for presentation.

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

## Manuscript methodology updates

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
