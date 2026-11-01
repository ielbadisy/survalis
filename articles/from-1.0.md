# Moving to the survalis verb interface

survalis 1.3.0 introduced a seven-verb interface. It keeps every prior
function working – nothing you have written breaks. The granular helpers
stay exported as the low-level layer under the verbs (see
`?survalis-lowlevel`). This vignette shows the verb that covers each of
them.

## The mapping

| Granular (still works) | Verb interface (recommended) |
|----|----|
| `fit_ranger(f, data, num.trees = 500)` | `fit(f, data, "ranger", spec = list(num.trees = 500))` |
| `predict_ranger(m, nd, times)` | `predict(m, nd, times)` |
| `predict_*()` then [`survmat_to_rmst()`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md) | `predict(m, nd, times, type = "rmst", tau = ...)` |
| `predict_*()` then [`survmat_to_quantile()`](https://elbadisyimad.com/survalis/reference/survmat_to_quantile.md) | `predict(m, nd, times, type = "median")` / `type = "quantile"` |
| `cv_survlearner(f, data, fit_ranger, predict_ranger, times, folds = 5)` | `evaluate(f, data, "ranger", times = ..., resampling = cv(5))` |
| `score_survmodel(m, times, metrics)` | `evaluate(m, times = ..., resampling = cv())` |
| `tune_ranger(f, data, times, param_grid, refit_best = TRUE)` | `tune(f, data, "ranger", grid = param_grid, times = ...)` |
| `benchmark(f, data, learners, times)` | `compare(f, data, models = learners, times = ...)` |
| `best_survlearner(bench, metric)` | [`print()`](https://rdrr.io/r/base/print.html) a `survalis_compare` (best is flagged) |
| `compute_pdp(m, data, "karno", times)` | `interpret(m, "pdp", feature = "karno", times = ...)` |
| `compute_shap(m, nd, data, times)` | `interpret(m, "shap", newdata = nd, times = ...)` |
| `compute_varimp(m, times)` | `interpret(m, "varimp", times = ...)` |
| `plot_pdp(res, "karno")` | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) a `survalis_interpret` |
| `estimate_rmst(m, data, tau, times, trt_col = "trt")` | `estimate(f, data, "<model>", treatment = "trt", tau = ...)` |

## Side-by-side

``` r

library(survalis)
f <- Surv(time, status) ~ age + karno + celltype

# 1.x
m1 <- fit_ranger(f, veteran, num.trees = 200)
p1 <- predict_ranger(m1, veteran[1:3, ], times = c(90, 180))

# 2.0
m2 <- fit(f, veteran, "ranger", spec = list(num.trees = 200), seed = 1)
p2 <- predict(m2, veteran[1:3, ], times = c(90, 180))

round(as.matrix(p2), 3)
#>       t=90 t=180
#> [1,] 0.404 0.269
#> [2,] 0.857 0.484
#> [3,] 0.649 0.435
```

## What changed on purpose

- **[`predict()`](https://rdrr.io/r/stats/predict.html) always returns a
  `survmat`** for the curve types: an `n x length(times)` matrix with
  the grid in `attr(x, "times")`. The `"t=<time>"` column names are
  still there, so the `*_survmat()` metric helpers keep accepting it.
- **`times` is optional.** When omitted,
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md)
  builds a grid over the observed event times and records it on the
  result (`res$times`), so a run is reproducible.
- **[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  /
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  take a resampling object**
  ([`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md))
  instead of a bare `folds` integer, and
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  evaluates every model on the *same* folds.
- **[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) engine
  arguments go in `spec`,** not `...`.
