# Hyperparameter Tuning for ranger Survival Models

Performs grid search tuning of a survival random forest model using
ranger over a set of hyperparameter combinations.

## Usage

``` r
tune_ranger(
  formula,
  data,
  times,
  param_grid = expand.grid(num.trees = c(100, 300), mtry = c(1, 2, 3), min.node.size =
    c(3, 5)),
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

  A `data.frame` containing the variables in the model.

- times:

  Numeric vector of time points at which to evaluate performance.

- param_grid:

  A `data.frame` of hyperparameter combinations. Must contain columns
  `num.trees`, `mtry`, and `min.node.size`.

- metrics:

  Character vector of evaluation metrics (e.g., `"cindex"`, `"ibs"`).

- folds:

  Number of cross-validation folds.

- seed:

  Random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the model using the best parameters.

- ...:

  Additional arguments passed to
  [`fit_ranger`](https://elbadisyimad.com/survalis/reference/fit_ranger.md).

## Value

If `refit_best = FALSE`, returns a `data.frame` of tuning results sorted
by the primary metric. If `refit_best = TRUE`, returns an
`"mlsurv_model"` object with the best parameters and an attribute
`"tuning_results"`.

## Details

Uses
[`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
to perform cross-validation for each parameter combination. If
`refit_best = TRUE`, the function returns the best-fitting model;
otherwise, it returns a tuning results table.

## See also

[`fit_ranger`](https://elbadisyimad.com/survalis/reference/fit_ranger.md),
[`predict_ranger`](https://elbadisyimad.com/survalis/reference/predict_ranger.md)

## Examples

``` r
mod_ranger_best <- tune_ranger(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(80, 160),
  param_grid = expand.grid(
    num.trees = c(25),
    mtry = c(2),
    min.node.size = c(5)
  ),
  metrics = c("cindex", "ibs"),
  folds = 2,
  refit_best = TRUE
)
summary(mod_ranger_best)
#> 
#> ── ranger summary ──────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: ranger
#> Learner: ranger
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
