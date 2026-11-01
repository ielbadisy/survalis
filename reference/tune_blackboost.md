# Tune blackboost Hyperparameters (Cross-Validation)

Cross-validates mboost Cox boosting models over a hyperparameter grid
and selects the best configuration according to the primary metric.
Optionally refits the best model on the full dataset.

## Usage

``` r
tune_blackboost(
  formula,
  data,
  times,
  param_grid = expand.grid(mstop = c(50, 100, 200), nu = c(0.05, 0.1), maxdepth = c(2,
    3)),
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

  A `data.frame` containing variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points (same scale as the training
  time).

- param_grid:

  A `data.frame` (e.g., from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) with
  columns such as `mstop`, `nu`, and `maxdepth`.

- metrics:

  Character vector of metrics to compute (e.g., `"cindex"`, `"ibs"`).
  The first entry is the primary selection metric.

- folds:

  Integer; number of CV folds. Default `5`.

- seed:

  Integer random seed for reproducibility. Default `123`.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns it. Default `FALSE`.

- ...:

  Additional arguments forwarded to
  [`fit_blackboost()`](https://elbadisyimad.com/survalis/reference/fit_blackboost.md).

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) of grid
results with metric columns, sorted by the primary metric. If
`refit_best = TRUE`, a fitted `mlsurv_model` from
[`fit_blackboost()`](https://elbadisyimad.com/survalis/reference/fit_blackboost.md)
using the selected hyperparameters.

## Details

Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_blackboost()`](https://elbadisyimad.com/survalis/reference/fit_blackboost.md)/[`predict_blackboost()`](https://elbadisyimad.com/survalis/reference/predict_blackboost.md)
so tuning mirrors the production path. Typical grids vary `mstop`, `nu`,
and tree `maxdepth`.

## Examples

``` r
# \donttest{
res_blackboost <- tune_blackboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(5, 10, 40),
  param_grid = expand.grid(
    mstop = c(100, 200),
    nu = c(0.05, 0.1),
    maxdepth = c(2, 3)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3
)
print(res_blackboost)
#> # basetable: 8 x 5
#>   mstop   nu maxdepth cindex   ibs
#> 1   100 0.05        2  0.713 0.072
#> 2   200 0.05        2  0.713 0.072
#> 3   100 0.10        2  0.713 0.072
#> 4   200 0.10        2  0.713 0.072
#> 5   100 0.05        3  0.713 0.072
#> 6   200 0.05        3  0.713 0.072
#> 7   100 0.10        3  0.713 0.072
#> 8   200 0.10        3  0.713 0.072

mod_blackboost_best <- tune_blackboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(5, 10, 40),
  param_grid = expand.grid(
    mstop = c(100, 200),
    nu = c(0.05, 0.1),
    maxdepth = c(2, 3)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)
summary(mod_blackboost_best)
#> 
#> ── blackboost summary ──────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: mboost
#> Learner: blackboost
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
# }
```
