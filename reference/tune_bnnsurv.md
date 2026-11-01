# Tune Bagged Nearest-Neighbour Survival Hyperparameters (Cross-Validation)

Cross-validates
[`fit_bnnsurv`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md)
over a user-supplied grid and aggregates metrics (for example,
`"cindex"`, `"ibs"`). Optionally refits and returns the best
configuration.

## Usage

``` r
tune_bnnsurv(
  formula,
  data,
  times,
  param_grid = expand.grid(k = c(2, 3), num_base_learners = c(30, 50), sample_fraction =
    c(0.5, 1), stringsAsFactors = FALSE),
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

  Training data frame.

- times:

  Numeric vector of evaluation time points used during CV.

- param_grid:

  A data frame (for example from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)) with
  columns `k`, `num_base_learners`, and `sample_fraction`. Defaults
  cover a small illustrative search.

- metrics:

  Character vector of metric names to compute and summarize. The first
  entry is treated as the primary ranking metric.

- folds:

  Integer number of cross-validation folds. Default `5`.

- seed:

  Integer random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns that model. Otherwise returns the results table.

## Value

If `refit_best = FALSE`, a `tibble` with one row per configuration
containing `metrics` and the tuning parameters, with a `failed` flag for
combinations that errored during CV. If `refit_best = TRUE`, an
`mlsurv_model` refit at the selected hyperparameters.

## Details

Evaluation is performed by
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
using
[`fit_bnnsurv()`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md)
and
[`predict_bnnsurv()`](https://elbadisyimad.com/survalis/reference/predict_bnnsurv.md)
to match production code paths. Results are ordered by the first entry
in `metrics`.

## See also

[`fit_bnnsurv()`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md),
[`predict_bnnsurv()`](https://elbadisyimad.com/survalis/reference/predict_bnnsurv.md),
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)

## Examples

``` r
# \donttest{
grid <- expand.grid(
  k = c(2, 3),
  num_base_learners = c(5, 10),
  sample_fraction = c(0.5, 1),
  stringsAsFactors = FALSE
)

res <- tune_bnnsurv(
  formula = Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran,
  times = c(100, 200),
  param_grid = grid,
  refit_best = FALSE
)
res
#> # basetable: 8 x 6
#>   k num_base_learners sample_fraction failed cindex   ibs
#> 1 2                 5             0.5  FALSE  0.567 0.271
#> 2 2                10             0.5  FALSE  0.567 0.251
#> 3 3                10             0.5  FALSE  0.566 0.251
#> 4 3                 5             0.5  FALSE  0.551 0.269
#> 5 2                 5             1.0  FALSE  0.529 0.300
#> 6 2                10             1.0  FALSE  0.529 0.287
#> 7 3                 5             1.0  FALSE  0.527 0.287
#> 8 3                10             1.0  FALSE  0.519 0.279

best <- tune_bnnsurv(
  formula = Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran,
  times = c(100, 200),
  param_grid = grid,
  refit_best = TRUE
)
head(predict_bnnsurv(best, newdata = veteran[1:5, ], times = c(50, 100)))
#>   t=50 t=100
#> 1  1.0   0.4
#> 2  0.7   0.5
#> 3  0.8   0.4
#> 4  0.8   0.8
#> 5  0.9   0.9
# }
```
