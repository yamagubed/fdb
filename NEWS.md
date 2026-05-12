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
