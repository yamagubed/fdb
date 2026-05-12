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
| `P2_GatedL1`        | \( \lambda |\delta| \cdot g(|\delta|/\widehat{SE}) \)| Smooth evidence-gated L1                        |
| `P3_SEScaledMCP`    | \( MCP(\delta; \lambda/\widehat{SE}, \gamma_{MCP}) \)| Information-adaptive minimax concave penalty    |
| `P4_LRWeightedL1`   | \( \lambda |\delta| \exp\{-\tfrac{1}{2}LR_0\} \)     | Likelihood-ratio-weighted L1                    |

All penalties operate on a smoothed absolute-value primitive
\(|\delta|_\varepsilon = \sqrt{\delta^2 + \varepsilon^2}\), so the resulting
penalized objective is twice differentiable in \(\delta\). This makes plug-in
sandwich variance estimation straightforward.

## Installation

From a local source build:

```r
# install.packages("remotes")
remotes::install_local("fdb_0.1.0.tar.gz", dependencies = TRUE)
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
control type I error under a prespecified set of population drift scenarios:

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
