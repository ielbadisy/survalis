# Tune a Survival Tree Model (`rpart`) via Cross-Validation

Performs hyperparameter tuning for the `rpart` survival tree model using
cross-validation, returning either a table of results or the best fitted
model.

## Usage

``` r
tune_rpart(
  formula,
  data,
  times,
  param_grid = expand.grid(minsplit = c(10, 20), cp = c(0.001, 0.01), maxdepth = c(10,
    30)),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1,
  refit_best = TRUE,
  ...
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the variables in the model.

- times:

  Numeric vector of time points for evaluation.

- param_grid:

  A `data.frame` or
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) result
  specifying hyperparameter combinations. Must include `minsplit`, `cp`,
  and `maxdepth`.

- metrics:

  Character vector of evaluation metrics to compute (e.g., `"cindex"`,
  `"ibs"`).

- folds:

  Number of cross-validation folds. Default is `5`.

- seed:

  Random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, the best hyperparameter combination is refitted on
  the full dataset.

- ...:

  Additional arguments passed to
  [`fit_rpart()`](https://elbadisyimad.com/survalis/reference/fit_rpart.md).

## Value

If `refit_best = FALSE`, returns a tibble summarizing mean CV
performance for each parameter combination. If `refit_best = TRUE`,
returns an `"mlsurv_model"` object for the best parameters, with a
`"tuning_results"` attribute.

## Examples

``` r

res_rpart <- tune_rpart(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    minsplit = c(10, 20),
    cp = c(0.001, 0.01),
    maxdepth = c(10, 30)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  seed = 42,
  refit_best = FALSE
)

print(res_rpart)
#> # basetable: 8 x 5
#>   minsplit    cp maxdepth cindex   ibs
#> 1       20 0.001       10  0.709 0.182
#> 2       20 0.001       30  0.709 0.182
#> 3       10 0.010       10  0.685 0.184
#> 4       20 0.010       10  0.685 0.184
#> 5       10 0.010       30  0.685 0.184
#> 6       20 0.010       30  0.685 0.184
#> 7       10 0.001       10  0.672 0.183
#> 8       10 0.001       30  0.672 0.183

mod_rpart_best <- tune_rpart(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    minsplit = c(10, 20),
    cp = c(0.001, 0.01),
    maxdepth = c(10, 30)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  seed = 42,
  refit_best = TRUE
)
```
