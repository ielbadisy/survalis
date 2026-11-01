# Tune XGBoost Survival Hyperparameters (Cross-Validation)

Cross-validates XGBoost survival models over a user-specified grid and
returns a results table with metric summaries per configuration. Any row
that errors during CV is marked `failed = TRUE`.

## Usage

``` r
tune_xgboost(
  formula,
  data,
  times,
  param_grid = expand.grid(nrounds = c(50, 100), max_depth = c(3, 6), eta = c(0.01, 0.1),
    aft_loss_distribution = c("extreme", "logistic"), aft_loss_distribution_scale =
    c(0.5, 1), objective = "survival:aft", stringsAsFactors = FALSE),
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
  columns: `nrounds`, `max_depth`, `eta`, `aft_loss_distribution`,
  `aft_loss_distribution_scale`, and `objective`. Values are passed
  through to
  [`fit_xgboost`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md)
  and ultimately to xgboost.

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

A `tibble` with one row per grid configuration, containing:

- nrounds, max_depth, eta, aft_loss_distribution,
  aft_loss_distribution_scale, objective:

  The grid values.

- failed:

  Logical; `TRUE` if the configuration errored.

- metric columns:

  One column per entry in `metrics` (when available).

The table is arranged by the first metric in `metrics` (ascending, as
implemented).

## Details

Internally calls
[`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_xgboost`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md)
/
[`predict_xgboost`](https://elbadisyimad.com/survalis/reference/predict_xgboost.md).
Any configuration that errors (e.g., due to invalid parameters or data
issues) is recorded with `failed = TRUE` and omitted from metric
summarization.

## See also

[`fit_xgboost`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md),
[`predict_xgboost`](https://elbadisyimad.com/survalis/reference/predict_xgboost.md)

## Examples

``` r
# \donttest{
  grid <- expand.grid(
    nrounds = c(20, 40),
    max_depth = c(2, 3),
    eta = c(0.1, 0.2),
    aft_loss_distribution = c("extreme", "logistic"),
    aft_loss_distribution_scale = c(0.5, 1),
    objective = "survival:aft",
    stringsAsFactors = FALSE
  )

  res_xgb <- tune_xgboost(
    formula = survival::Surv(time, status) ~ age + karno + celltype,
    data    = survival::veteran,
    times   = c(100, 200),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds   = 2,
    seed    = 123
  )
  head(res_xgb)
#> # basetable: 6 x 9
#>   nrounds max_depth eta aft_loss_distribution aft_loss_distribution_scale
#> 1      20         2 0.1              logistic                         0.5
#> 2      40         2 0.1              logistic                         0.5
#> 3      40         2 0.2               extreme                         0.5
#> 4      40         3 0.1              logistic                         0.5
#> 5      20         3 0.1              logistic                         0.5
#> 6      40         2 0.2               extreme                         1.0
#>      objective failed cindex   ibs
#> 1 survival:aft  FALSE  0.706 0.186
#> 2 survival:aft  FALSE  0.696 0.184
#> 3 survival:aft  FALSE  0.690 0.198
#> 4 survival:aft  FALSE  0.678 0.188
#> 5 survival:aft  FALSE  0.677 0.186
#> 6 survival:aft  FALSE  0.676 0.200
# }
```
