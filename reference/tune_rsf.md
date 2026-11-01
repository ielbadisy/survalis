# Tune Random Survival Forest Hyperparameters (Cross-Validation)

Cross-validates RSF models over a specified parameter grid and selects
the best configuration according to the primary metric. Optionally
refits the best model on the full dataset.

## Usage

``` r
tune_rsf(
  formula,
  data,
  times,
  param_grid = expand.grid(ntree = c(200, 500), mtry = c(1, 2, 3), nodesize = c(5, 15)),
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

  Numeric vector of evaluation time points (same scale as the survival
  time used at training). Must be non-negative and finite.

- param_grid:

  A data frame (e.g., from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) containing
  columns for RSF hyperparameters such as `ntree`, `mtry`, and
  `nodesize`.

- metrics:

  Character vector of metrics to evaluate/optimize (e.g., `"cindex"`,
  `"ibs"`). The first entry is the primary selection metric.

- folds:

  Integer; number of cross-validation folds.

- seed:

  Integer random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns it, with tuning results attached.

- ...:

  Additional arguments forwarded to the underlying engine where
  applicable.

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) of grid
results with metric columns and hyperparameters, ordered by the first
metric. If `refit_best = TRUE`, a fitted `mlsurv_model` from
[`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md)
with attribute `"tuning_results"` containing the full grid results.

## Details

Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md)/[`predict_rsf()`](https://elbadisyimad.com/survalis/reference/predict_rsf.md)
so tuning uses the same code paths as production. Typical grids vary
`ntree`, `mtry`, and `nodesize`.

## See also

[`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md),
[`predict_rsf()`](https://elbadisyimad.com/survalis/reference/predict_rsf.md)

## Examples

``` r
# \donttest{
mod_rsf_best <- tune_rsf(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    ntree = c(200, 500),
    mtry = c(1, 2, 3),
    nodesize = c(5, 15)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

summary(mod_rsf_best)
#> 
#> ── rsf summary ─────────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + celltype + karno
#> Engine: randomForestSRC
#> Learner: rsf
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, celltypesmallcell, celltypeadeno, celltypelarge, karno"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
# }
```
