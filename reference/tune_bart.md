# Tune BART Survival Hyperparameters (Cross-Validation)

Cross-validates BART survival models over a user‐supplied grid and
selects the best configuration by the primary metric. Optionally refits
the best model on the full dataset.

## Usage

``` r
tune_bart(
  formula,
  data,
  times,
  param_grid = expand.grid(K = c(3, 5), ntree = c(50, 100), power = c(2, 2.5), base =
    c(0.75, 0.95)),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1,
  refit_best = FALSE
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` with variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points.

- param_grid:

  A `data.frame` (e.g., from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) of
  candidate hyperparameters (e.g., `K`, `ntree`, `power`, `base`).

- metrics:

  Character vector of evaluation metrics (default `c("cindex","ibs")`).
  The first entry is used as the primary selection metric.

- folds:

  Integer; number of cross‐validation folds (default `5`).

- seed:

  Integer random seed for reproducibility (default `123`).

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns it. If `FALSE`, returns a results table.

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) sorted by
the primary metric with one row per grid combination. If
`refit_best = TRUE`, a fitted `mlsurv_model` returned by
[`fit_bart`](https://elbadisyimad.com/survalis/reference/fit_bart.md).

## Details

Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_bart()`](https://elbadisyimad.com/survalis/reference/fit_bart.md)/[`predict_bart()`](https://elbadisyimad.com/survalis/reference/predict_bart.md)
so tuning mirrors the production prediction path.

## See also

[`fit_bart`](https://elbadisyimad.com/survalis/reference/fit_bart.md),
[`predict_bart`](https://elbadisyimad.com/survalis/reference/predict_bart.md),
[`surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html)

## Examples

``` r
# \donttest{
grid <- expand.grid(K = c(3), ntree = c(50), power = c(2), base = c(0.75, 0.95))
res <- tune_bart(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  param_grid = grid,
  times = c(10, 60),
  refit_best = FALSE
)
print(res)
#> # basetable: 2 x 6
#>   K ntree power base cindex   ibs
#> 1 3    50     2 0.75  0.696 0.152
#> 2 3    50     2 0.95  0.696 0.152

mod_bart <- tune_bart(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  param_grid = grid,
  times = c(10, 60),
  refit_best = TRUE
)
summary(mod_bart)
#> 
#> ── bart summary ────────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: BART
#> Learner: bart
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
# }
```
