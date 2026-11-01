# Tune CoxBoost Survival Hyperparameters (Cross-Validation)

Cross-validates CoxBoost survival models over a user-specified grid and
returns a results table with metric summaries per configuration. Any row
that errors during CV is marked `failed = TRUE`.

## Usage

``` r
tune_coxboost(
  formula,
  data,
  times,
  param_grid = expand.grid(stepno = c(50, 100), penalty = c(100, 200), criterion =
    "pscore", stepsize.factor = 1, stringsAsFactors = FALSE),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points.

- param_grid:

  A `data.frame` (e.g., from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) with
  columns: `stepno`, `penalty`, `criterion`, `stepsize.factor`. Values
  are passed through to
  [`fit_coxboost`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md).

- metrics:

  Character vector of metrics to evaluate (e.g., `"cindex"`, `"ibs"`).
  The first entry is treated as the primary selection metric for
  ordering.

- folds:

  Integer number of CV folds.

- seed:

  Integer random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

## Value

A `data.frame` with one row per grid configuration, containing the grid
values, a `failed` column, and one column per metric.

## Details

Internally calls
[`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_coxboost`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md)
/
[`predict_coxboost`](https://elbadisyimad.com/survalis/reference/predict_coxboost.md).
Any configuration that errors is recorded with `failed = TRUE` and
omitted from metric summarization.

## See also

[`fit_coxboost`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md),
[`predict_coxboost`](https://elbadisyimad.com/survalis/reference/predict_coxboost.md)

## Examples

``` r
# \donttest{
  grid <- expand.grid(
    stepno = c(50, 100),
    penalty = c(100, 200),
    criterion = "pscore",
    stepsize.factor = 1,
    stringsAsFactors = FALSE
  )

  res_cb <- tune_coxboost(
    formula = survival::Surv(time, status) ~ age + karno + celltype,
    data    = survival::veteran,
    times   = c(100, 200),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds   = 2,
    seed    = 123
  )
  head(res_cb)
#> # basetable: 4 x 7
#>   stepno penalty criterion stepsize.factor failed cindex   ibs
#> 1     50     200    pscore               1  FALSE  0.730 0.162
#> 2    100     200    pscore               1  FALSE  0.730 0.161
#> 3     50     100    pscore               1  FALSE  0.729 0.162
#> 4    100     100    pscore               1  FALSE  0.729 0.162
# }
```
