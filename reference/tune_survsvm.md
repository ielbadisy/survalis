# Tune Survival SVM Hyperparameters (Cross-Validation)

Cross-validates Survival SVM models over a user-specified grid and
selects the best configuration based on the chosen metric. Optionally
refits the best model on the full dataset.

## Usage

``` r
tune_survsvm(
  formula,
  data,
  times,
  metrics = "cindex",
  param_grid = list(gamma.mu = c(0.01, 0.1, 1), kernel = c("lin_kernel", "add_kernel")),
  folds = 5,
  seed = 42,
  ncores = 1,
  refit_best = FALSE,
  dist = "km",
  shape = 1
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`. The
  left-hand side must be a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object from the
  `survival` package.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points (same scale as the survival
  time used at training). Must be non-negative and finite.

- metrics:

  Character vector of metrics to evaluate/optimize (e.g., `"cindex"`,
  `"ibs"`). The first entry is the primary selection metric.

- param_grid:

  A named list of candidate hyperparameters with entries `gamma.mu` and
  `kernel`.

- folds:

  Number of cross-validation folds.

- seed:

  Integer random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best configuration on the full data and
  returns it as the result (with tuning attributes).

- dist, shape:

  Passed to
  [`predict_survsvm()`](https://elbadisyimad.com/survalis/reference/predict_survsvm.md)
  during cross-validation.

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tune_surv"`) of grid
results with metric columns and tuning parameters. If
`refit_best = TRUE`, a fitted `mlsurv_model` (class augmented with
`"tune_surv"`) with the full results attached in
`attr(, "tuning_results")`.

## Details

Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_survsvm()`](https://elbadisyimad.com/survalis/reference/fit_survsvm.md)/[`predict_survsvm()`](https://elbadisyimad.com/survalis/reference/predict_survsvm.md)
to keep code paths consistent with production usage. A grid point that
fails is skipped with a warning.

## See also

[`fit_survsvm()`](https://elbadisyimad.com/survalis/reference/fit_survsvm.md),
[`predict_survsvm()`](https://elbadisyimad.com/survalis/reference/predict_survsvm.md)

## Examples

``` r
# \donttest{
grid <- list(
  gamma.mu = c(0.01, 0.1),
  kernel = c("lin_kernel", "add_kernel")
)

res_svm <- tune_survsvm(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs"),
  param_grid = grid,
  folds = 3
)
res_svm
#> # basetable: 4 x 4
#>   gamma.mu     kernel cindex   ibs
#> 1     0.10 lin_kernel  0.727 0.141
#> 2     0.10 add_kernel  0.720 0.148
#> 3     0.01 lin_kernel  0.681 0.159
#> 4     0.01 add_kernel  0.536 0.171
# }
```
