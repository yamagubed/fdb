# fdb

**Frequentist Dynamic Borrowing for Hybrid-Control Survival Trials**

The `fdb` package implements a class of likelihood-informed frequentist dynamic
borrowing methods for hybrid-control survival trials, based on penalized Cox
partial likelihood estimation. The borrowing strength is controlled by a single
tuning parameter \(\lambda\) acting on a drift parameter that captures
population shift between concurrent and external controls. Both standard
model-based standard errors and penalized estimating-equation sandwich standard
errors are provided.

## Penalty methods

The package implements four likelihood-informed penalties together with the
adaptive lasso comparator of Li et al. (2023):

| Tag                 | Penalty                                              | Description                                     |
| ------------------- | ---------------------------------------------------- | ----------------------------------------------- |
| `LiAdaptiveLasso`   | \( \lambda |\hat\delta_0|^{-\gamma} |\delta| \)      | Adaptive lasso (Li et al. 2023)                 |
| `P1_SEScaledL1`     | \( \lambda |\delta| / \widehat{SE}(\hat\delta_0) \)  | Precision-weighted L1                           |
| `P2_GatedL1`        | \( \lambda\int_0^{|\delta|}g(u/\widehat{SE})\,du \)| Smoothed integrated-gate                        |
| `P3_SEScaledMCP`    | \( MCP(\delta; \lambda/\widehat{SE}, \gamma_{MCP}) \)| Information-adaptive minimax concave penalty    |
| `P4_LRWeightedL1`   | \( \lambda |\delta| \exp\{-\tfrac{1}{2}LR_0\} \)     | Likelihood-ratio-weighted L1                    |

All penalties operate on a smoothed absolute-value primitive
\(|\delta|_\varepsilon = \sqrt{\delta^2 + \varepsilon^2}\), so the resulting
penalized objective is twice differentiable in \(\delta\). This makes plug-in
sandwich variance estimation possible as a local plug-in approximation.
First-stage estimates are held fixed; conditional or unconditional validity
is not guaranteed. Negative P2/P3 curvature is clipped only for the variance
calculation; raw curvature remains available as `pen_curv_raw`.

P2 integrates the gate from `eps` to `sqrt(delta^2 + eps^2)`. MCP integrates
a C1 slope from zero to `sqrt(delta^2 + eps^2) - eps`, smoothing both the
origin and the flat-tail transition. Set `rho_mcp` in direct fits, calibration,
or `run_fdb_study`, or in the tuning list passed to `run_simulation`.
Its software default is 0.1; it is not a calibrated manuscript value.
Existing function names and method tags are retained for compatibility.
Both P2 and P3 results and calibrated lambdas from the old definitions
need recomputation. No-covariate studies are supported with `p = 0`.

## Installation

From a local source build:

```r
# install.packages("remotes")
remotes::install_local("fdb_0.2.0.tar.gz", dependencies = TRUE)
```

## Quick start

```r
library(fdb)

# Simulate a hybrid-control trial with HR = 0.8 and no population drift
set.seed(1)
sim <- simulate_hybrid_cox(
  nI1 = 150, nI0 = 150, nE = 300,
  theta0 = log(0.8), delta0 = 0
)

# Fit all methods on this dataset
fit_all_methods(sim$data)

# Fit a single method
fit_one_penalized_method(sim$data, method = "P1", lambda = 0.2)
```

## Design-stage calibration and operating characteristics

For confirmatory use, the penalty strength \(\lambda\) should be calibrated to
target an error-rate threshold over a prespecified drift grid. Finite Monte
Carlo calibration does not guarantee control between grid points or outside
that grid:

```r
study <- run_fdb_study(
  scenario_base   = scenario_S1,
  drift_hr_range  = c(0.8, 1.2),
  drift_by_hr     = 0.05,
  alpha           = 0.025,
  alt_hr          = 0.8,
  do_calibration  = TRUE,
  nsim_cal        = 500,
  nsim_curve      = 1000,
  parallel        = TRUE,
  ncores          = 2,
  seed            = 1
)

study$type1_curve
study$power_curve
study$calibration$lambda_star
```

## References

- Li, R., Lin, R., Huang, J., Tian, L., and Zhu, J. (2023).
  A frequentist approach to dynamic borrowing.
  *Biometrical Journal* 65(7), 2100406.
- Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
  concave penalty. *The Annals of Statistics* 38(2), 894-942.
- Andersen, P. K. and Gill, R. D. (1982). Cox's regression model for
  counting processes: A large sample study.
  *The Annals of Statistics* 10(4), 1100-1120.

## License

MIT (see the `LICENSE` file).

## Interpretation and reproducibility

The model-based SE conditions on the estimated drift as a fixed offset and
can substantially underestimate uncertainty. The sandwich SE is a local
plug-in approximation holding adaptive weights fixed; nominal coverage is
not guaranteed. Report observed coverage, bias, RMSE and valid-fit counts.
An `alpha_cal` threshold above `alpha` permits error inflation relative to
the nominal test level. ESS is a variance-equivalent gain, not a literal
number of borrowed controls; negative values indicate a precision loss and
positive values do not establish low bias.

Calibration never substitutes an uncalibrated default when it fails. Refine
the prespecified design or candidate grid and repeat confirmation. An explicit
`drift_set_confirm` is honored by the calibration functions; when omitted,
confirmation uses the calibration grid. The one-stop study wrapper confirms
on its calibration grid and evaluates performance on its wider curve grid.

`run_simulation()` always returns replicate estimates in `$raw`. For curves:

```r
# Larger nsim is needed for scientific conclusions.
curve <- run_drift_curve(theta0 = 0, drift_set = log(c(1, 1.1)),
                         scenario_base = scenario_S1, lambdas = lambdas_default,
                         nsim = 2, keep_raw = TRUE)
curve$summary
head(curve$raw)
```

The default `keep_raw = FALSE` preserves the summary-only curve interface.
`run_fdb_study(keep_raw = TRUE)` also returns `raw_type1` and `raw_power`.
Save these with `saveRDS()` for paired ESS uncertainty and normality diagnostics.
Simulation summaries include `n_valid`, `n_missing`, and Monte Carlo standard
errors. Check those counts before interpreting rates. Parallel execution is
opt-in and defaults to two workers. Install the package before starting
workers; fixed seeds are reproducible for a fixed worker count and environment,
not necessarily between serial runs and different parallel configurations.
