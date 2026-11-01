# Tune Oblique Random Survival Forests (ORSF) via Cross-Validation

Performs grid-search hyperparameter tuning for aorsf models using
cross-validation and selects the best configuration by the primary
metric. Optionally refits the best model on the full dataset.

## Usage

``` r
tune_orsf(
  formula,
  data,
  times,
  param_grid = expand.grid(n_tree = c(100, 300), mtry = c(2, 3), min_events = c(5, 10)),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1,
  refit_best = FALSE,
  ...
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points (same scale as the training
  survival time). Must be positive and finite.

- param_grid:

  A `data.frame` (e.g., from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) of
  candidate hyperparameters. Typical columns include:

  - `n_tree`: number of trees (e.g., `c(100, 300)`)

  - `mtry`: number of variables tried at each split

  - `min_events`: minimum number of events required to attempt a split

- metrics:

  Character vector of metrics to evaluate/optimize (e.g., `"cindex"`,
  `"ibs"`). The first entry is used as the primary selection metric.

- folds:

  Integer; number of cross-validation folds. Default is `5`.

- seed:

  Integer random seed for reproducibility. Default is `123`.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns it. Default is `FALSE`.

- ...:

  Additional arguments forwarded to the underlying engine via
  [`fit_orsf`](https://elbadisyimad.com/survalis/reference/fit_orsf.md).

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) with one
row per grid combination and columns for hyperparameters and metrics,
ordered by the first metric. If `refit_best = TRUE`, a fitted
`mlsurv_model` from
[`fit_orsf`](https://elbadisyimad.com/survalis/reference/fit_orsf.md)
using the best settings.

## Details

Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_orsf`](https://elbadisyimad.com/survalis/reference/fit_orsf.md) /
[`predict_orsf`](https://elbadisyimad.com/survalis/reference/predict_orsf.md)
so tuning uses the exact same code paths as production. The `min_events`
column of `param_grid` is passed to the engine as `n_split` (minimum
events for a candidate split).

## See also

[`fit_orsf`](https://elbadisyimad.com/survalis/reference/fit_orsf.md),
[`predict_orsf`](https://elbadisyimad.com/survalis/reference/predict_orsf.md),
aorsf

## Examples

``` r
# \donttest{
  res_orsf <- tune_orsf(
    formula    = Surv(time, status) ~ age + karno,
    data       = veteran,
    times      = c(100, 200, 300),
    param_grid = expand.grid(
      n_tree     = c(100, 200),
      mtry       = c(1, 2),
      min_events = c(5, 10)
    ),
    metrics    = c("cindex", "ibs"),
    folds      = 2
  )
  print(res_orsf)
#> # basetable: 8 x 5
#>   n_tree mtry min_events cindex   ibs
#> 1    100    1          5  0.652 0.179
#> 2    200    1          5  0.651 0.178
#> 3    200    1         10  0.647 0.180
#> 4    100    1         10  0.645 0.180
#> 5    100    2         10  0.634 0.184
#> 6    200    2          5  0.629 0.182
#> 7    200    2         10  0.627 0.182
#> 8    100    2          5  0.624 0.185

  mod_orsf_best <- tune_orsf(
    formula    = Surv(time, status) ~ age + karno,
    data       = veteran,
    times      = c(100, 200, 300),
    param_grid = expand.grid(
      n_tree     = c(100, 200),
      mtry       = c(1, 2),
      min_events = c(5, 10)
    ),
    metrics    = c("cindex", "ibs"),
    folds      = 2,
    refit_best = TRUE
  )

  summary(mod_orsf_best)
#> 
#> ── orsf summary ────────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno
#> Engine: aorsf
#> Learner: orsf
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
# }
```
