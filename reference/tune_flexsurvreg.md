# Tune Parametric Survival Models with flexsurvreg

Performs hyperparameter tuning for
[`fit_flexsurvreg`](https://elbadisyimad.com/survalis/reference/fit_flexsurvreg.md)
over a grid of parametric distributions using cross-validation. Returns
either a summary table of performance metrics or a refitted best model.

## Usage

``` r
tune_flexsurvreg(
  formula,
  data,
  times,
  param_grid = c("weibull", "exponential", "lognormal"),
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

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)-based survival
  formula.

- data:

  A data frame containing the variables in the model.

- times:

  Numeric vector of time points at which to evaluate performance.

- param_grid:

  Character vector of distributions to test (default =
  `c("weibull", "exponential", "lognormal")`).

- metrics:

  Character vector of performance metrics to compute (default =
  `c("cindex", "ibs")`).

- folds:

  Integer, number of cross-validation folds (default = 5).

- seed:

  Random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best-performing model on the full
  dataset and returns it. If `FALSE`, returns the tuning results.

- ...:

  Additional arguments passed to
  [`fit_flexsurvreg`](https://elbadisyimad.com/survalis/reference/fit_flexsurvreg.md).

## Value

If `refit_best = FALSE`, returns a data frame of tuning results with one
row per distribution and columns for each metric. If
`refit_best = TRUE`, returns the best `mlsurv_model` object.

## Examples

``` r

res <- tune_flexsurvreg(Surv(time, status) ~ age + karno + celltype,
                        data = veteran,
                        param_grid = c("weibull", "exponential", "lognormal"),
                        times = c(100, 200, 300))
print(res)
#> # basetable: 3 x 3
#>          dist cindex   ibs
#> 1   lognormal  0.726 0.156
#> 2 exponential  0.723 0.153
#> 3     weibull  0.721 0.153

best_mod <- tune_flexsurvreg(Surv(time, status) ~ age + karno + celltype,
                             data = veteran,
                             param_grid = c("weibull", "exponential", "lognormal"),
                             times = c(100, 200, 300),
                             refit_best = TRUE)
summary(best_mod)
#> 
#> ── flexsurvreg summary ─────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: flexsurv
#> Learner: flexsurvreg
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
