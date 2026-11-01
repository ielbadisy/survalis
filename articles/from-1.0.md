# Moving from survalis 1.x to the seven verbs

survalis 2.0.0 exports seven verbs and nothing that duplicates them:
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
and
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md).
Around them sit the resampling specs
([`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md))
and three lookups that give the valid ids
([`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md),
[`list_metrics()`](https://elbadisyimad.com/survalis/reference/list_metrics.md),
[`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md)).

The 1.x per-learner functions (`fit_ranger()`, `predict_ranger()`,
`tune_ranger()`, …) and the helper families (`cv_*()`, `benchmark*()`,
`compute_*()`, `plot_*()`, `*_survmat()`) are no longer exported. Each
one maps to a verb argument: the learner goes in `model`, the metric in
`metrics`, the explanation in `method`, the output scale in `type`, and
every result has a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

## The mapping

### Fitting and prediction

| 1.x | 2.0 |
|----|----|
| `fit_ranger(f, data, num.trees = 500)` | `fit(f, data, "ranger", spec = list(num.trees = 500))` |
| `predict_ranger(m, nd, times)` | `predict(m, nd, times)` |
| `survmat_to_rmst(S, times, tau)` | `predict(m, nd, times, type = "rmst", tau = ...)` |
| `survmat_to_quantile(S, times, p)` | `predict(m, nd, times, type = "quantile", probs = p)` or `type = "median"` |
| `survmat_to_chf(S)` | `predict(m, nd, times, type = "chf")` |
| `survmat_to_haz(S, times)` | `predict(m, nd, times, type = "hazard")` |
| `fit_survmetalearner(...)` | `fit(f, data, "survmetalearner", spec = list(learners = ...))` |
| a hand-written `fit_fun` / `pred_fun` pair | `fit(f, data, model = list(fit = , predict = ))` |

### Evaluation

| 1.x | 2.0 |
|----|----|
| `cv_survlearner(f, data, fit_fun, pred_fun, times, folds = 5)` | `evaluate(f, data, "ranger", times = ..., resampling = cv(5))` |
| `score_survmodel(m, times, metrics)` | `evaluate(m, times = ..., metrics = ...)` |
| `cv_summary(res)` | `summary(evaluate(...))` |
| `cv_plot(res)` | `plot(evaluate(...))` |
| `cindex_survmat(y, S, t)`, `ibs_survmat(y, S, times)`, `brier(...)`, `auc_survmat()`, `iae_survmat()`, `ise_survmat()`, `ece_survmat()` | `evaluate(S, y = y, metrics = c("cindex", "ibs", ...))` |
| predict on a test set, then a metric helper | `evaluate(m, newdata = test, times = ...)` |
| `roc_survmat()` + `plot_roc()` | `plot(interpret(m, "roc", newdata = test, eval_time = t))` |
| `dca_survmat()` + `plot_dca()` | `plot(interpret(m, "dca", newdata = test, eval_time = t))` |

### Tuning and comparison

| 1.x | 2.0 |
|----|----|
| `tune_ranger(f, data, times, param_grid, refit_best = TRUE)` | `tune(f, data, "ranger", grid = param_grid, times = ...)` |
| `list_tunable_survlearners()` | `list_survlearners(has_tune = TRUE)` |
| `benchmark(f, data, learners, times)` | `compare(f, data, models = learners, times = ...)` |
| `benchmark_default_survlearners(...)` | `compare(...)` |
| `benchmark_tuned_survlearners(...)` | `compare(..., tune = TRUE)` |
| `summarise_benchmark()`, `best_survlearner()` | [`print()`](https://rdrr.io/r/base/print.html) or [`summary()`](https://rdrr.io/r/base/summary.html) a [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md) result |
| `plot_benchmark()` | `plot(compare(...))` |

### Interpretation and causal contrasts

| 1.x | 2.0 |
|----|----|
| `compute_pdp(m, data, "karno", times)` + `plot_pdp()` | `plot(interpret(m, "pdp", feature = "karno", times = ...))` |
| `compute_ale()` + `plot_ale()` | `interpret(m, "ale", feature = ...)` |
| `compute_shap()` + `plot_shap()` | `interpret(m, "shap", newdata = nd, ...)` |
| `compute_shap_matrix()` + `plot_shap_beeswarm()` | `interpret(m, "shap_matrix", newdata = nd, ...)` |
| `compute_varimp()` + `plot_varimp()` | `interpret(m, "varimp", ...)` |
| `compute_interactions()` + `plot_interactions()` | `interpret(m, "interaction", ...)` |
| `compute_counterfactual()` | `interpret(m, "counterfactual", target_time = ...)` |
| `compute_surrogate()`, `compute_tree_surrogate()` | `interpret(m, "surrogate", ...)`, `interpret(m, "tree_surrogate", ...)` |
| `compute_calibration()` + `plot_calibration()` | `interpret(m, "calibration", eval_time = ...)` |
| `compute_survcluster()` + `plot_survcluster()` | `interpret(m, "survcluster", K = ...)` |
| `estimate_rmst(m, data, tau, times, trt_col = "trt")` + `plot_rmst()` | `estimate(f, data, "<model>", treatment = "trt", tau = ...)` |

### Removed without a replacement

- `screen_fdr()`: use the `"selectcox"` or `"glmnet"` learners, which
  select variables while fitting.
- `plot_survcurve()`, `theme_survalis()`, `scale_color_survalis()` and
  `scale_fill_survalis()`: survalis plots still use this style
  internally; for Kaplan-Meier figures use a dedicated package such as
  ‘ggsurvfit’.

## Side by side

``` r

library(survalis)
f <- Surv(time, status) ~ age + karno + celltype

m <- fit(f, veteran, "ranger", spec = list(num.trees = 200), seed = 1)
round(as.matrix(predict(m, veteran[1:3, ], times = c(90, 180))), 3)
#>       t=90 t=180
#> [1,] 0.404 0.269
#> [2,] 0.857 0.484
#> [3,] 0.649 0.435
predict(m, veteran[1:3, ], times = c(90, 180), type = "rmst")
#>         1         2         3 
#>  93.42357 143.92197 122.97001
```

A learner of your own replaces the 1.x `fit_fun` / `pred_fun` arguments.
Its `predict` returns one row per subject and one column per time:

``` r

expo <- list(
  name = "exponential",
  fit = function(formula, data, ...) survival::survreg(formula, data, dist = "exponential"),
  predict = function(object, newdata, times) {
    rate <- 1 / predict(object, newdata, type = "response")
    exp(-outer(rate, times))
  }
)
evaluate(Surv(time, status) ~ age + karno, veteran, model = expo,
         times = c(90, 180), resampling = cv(3, seed = 1))
#> <survalis_eval> model = exponential
#> <survalis_resampling> 3-fold CV, strata = status, seed = 1
#>   grid: 2 times in [90, 180]
#>   cindex  0.7070  (sd 0.0420, 95% [0.6590, 0.7540], n = 3)
#>   ibs     0.1810  (sd 0.0120, 95% [0.1680, 0.1950], n = 3)
```

## Other changes on purpose

- **[`predict()`](https://rdrr.io/r/stats/predict.html) returns a
  survival matrix** for the curve types: an `n x length(times)` matrix
  with the grid in `attr(x, "times")` and column names `"t=<time>"`.
- **`times` is optional.** When omitted, a grid of up to 50 event-time
  quantiles (5th to 95th percentile) is used and recorded on the result
  (`res$times`), so a run is reproducible.
- **Resampling is an object**
  ([`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md))
  instead of a bare `folds` integer, and
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  evaluates every model on the *same* splits.
- **Engine arguments go in `spec`,** not `...`.
