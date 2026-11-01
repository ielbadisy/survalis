# Tune a survival learner over a hyperparameter grid

Grid search over resampling folds. Learners listed by
`list_survlearners(has_tune = TRUE)` have a built-in default grid; any
other learner, including a custom one (see
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md)), is tuned
over the `grid` you supply, with every configuration scored on the same
splits. The best configuration is refit on the full data and returned as
a ready-to-use `survalis_fit`.

## Usage

``` r
tune(
  formula,
  data,
  model,
  grid = NULL,
  times = NULL,
  resampling = cv(),
  metric = "ibs",
  ncores = 1,
  ...
)
```

## Arguments

- formula:

  A survival formula.

- data:

  A data frame.

- model:

  A learner id (`list_survlearners()$learner`) or a custom learner.

- grid:

  A data frame of hyperparameter combinations (one row per
  configuration, for example from
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)). `NULL`
  uses the learner's built-in default grid, if it has one.

- times:

  Evaluation grid; `NULL` picks one from the data quantiles.

- resampling:

  A resampling spec. Built-in grids use
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  only (its `v` and `seed` drive the inner cross-validation); a supplied
  `grid` for a learner without a built-in tuner accepts any spec.

- metric:

  A single metric to optimise (default `"ibs"`).

- ncores:

  Cores for the inner CV (default `1`).

- ...:

  Forwarded to the learner's built-in tuner.

## Value

A `survalis_tune` object: `fit` (the refit `survalis_fit` at the best
configuration), `results` (the scored grid), `best` (one-row best
config), `times`, `resampling`, `metric`, `model`.

## Examples

``` r
tuned <- tune(Surv(time, status) ~ age + karno, veteran, "ranger",
              grid = expand.grid(num.trees = c(100, 300), mtry = 2, min.node.size = c(5, 15)),
              times = c(90, 180), resampling = cv(3), metric = "ibs")
tuned
#> <survalis_tune> model = ranger  |  optimise ibs over 4 configs
#> <survalis_resampling> 3-fold CV, strata = status
#>   best: num.trees = 100, mtry = 2, min.node.size = 15
#>   ibs at best: 0.2010
predict(tuned, veteran[1:3, ], times = c(90, 180))
#> <survmat> survival  |  3 obs x 2 times  [90, 180]
#>           t=90     t=180
#> [1,] 0.3742667 0.1699102
#> [2,] 0.8498103 0.3458225
#> [3,] 0.5941357 0.3412420
```
