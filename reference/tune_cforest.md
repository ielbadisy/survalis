# Tune a Conditional Inference Survival Forest

Performs cross-validated hyperparameter tuning for a conditional
inference survival forest using the
[`fit_cforest()`](https://elbadisyimad.com/survalis/reference/fit_cforest.md)
and
[`predict_cforest()`](https://elbadisyimad.com/survalis/reference/predict_cforest.md)
functions.

## Usage

``` r
tune_cforest(
  formula,
  data,
  times,
  param_grid = expand.grid(ntree = c(100, 300), mtry = c(2, 3), mincriterion = c(0,
    0.95), fraction = c(0.5, 0.632)),
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

  A [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) survival
  formula specifying the outcome and predictors.

- data:

  A data frame containing the variables in the model.

- times:

  A numeric vector of time points at which to evaluate performance.

- param_grid:

  A data frame or list specifying the grid of hyperparameter values to
  evaluate. Columns should include `ntree`, `mtry`, `mincriterion`, and
  `fraction`.

- metrics:

  A character vector of performance metrics to compute (default =
  `c("cindex", "ibs")`).

- folds:

  Integer, number of cross-validation folds (default = `5`).

- seed:

  Integer, random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical, if `TRUE` refits a model with the best-performing
  hyperparameters and returns it; otherwise returns a summary table of
  results.

- ...:

  Additional arguments passed to
  [`fit_cforest()`](https://elbadisyimad.com/survalis/reference/fit_cforest.md).

## Value

If `refit_best = FALSE`, a data frame of mean cross-validation scores
for each hyperparameter combination (class `"tuned_surv"`). If
`refit_best = TRUE`, an `mlsurv_model` object fitted with the optimal
hyperparameters.

## Details

Cross-validation is performed using
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
and results are sorted in descending order of the first metric specified
in `metrics`.

## Examples

``` r
grid <- expand.grid(
  ntree = c(50, 100),
  mtry = c(2, 4),
  mincriterion = c(0, 0.95),
  fraction = c(0.632)
)
res <- tune_cforest(Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 200),
  param_grid = grid,
  folds = 3
)
print(res)
#> # basetable: 8 x 6
#>   ntree mtry mincriterion fraction cindex   ibs
#> 1    50    2         0.95    0.632  0.713 0.173
#> 2   100    4         0.95    0.632  0.713 0.161
#> 3   100    2         0.95    0.632  0.706 0.171
#> 4    50    4         0.95    0.632  0.704 0.166
#> 5    50    4         0.00    0.632  0.691 0.159
#> 6   100    4         0.00    0.632  0.690 0.158
#> 7   100    2         0.00    0.632  0.672 0.161
#> 8    50    2         0.00    0.632  0.668 0.167
```
