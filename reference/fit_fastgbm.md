# Fit a fastgbm Survival Model

Fits a histogram-based gradient boosting survival model with the fastgbm
package and returns an object compatible with the `mlsurv_model`
interface. The piecewise-exponential objective (`objective = "pexp"`) is
used by default because it models the hazard jointly over covariates and
time, so predicted survival curves can be evaluated at arbitrary time
points.

## Usage

``` r
fit_fastgbm(
  formula,
  data,
  objective = c("pexp", "cox", "aft"),
  ntrees = 200,
  learning_rate = 0.1,
  max_depth = 5,
  pexp_bins = 10,
  ...
)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the formula variables.

- objective:

  fastgbm survival objective: `"pexp"` (default), `"cox"`, or `"aft"`.
  Only `"pexp"` supports `predict(type = "survival")` at arbitrary times
  and is required for the survalis prediction contract.

- ntrees:

  Number of boosting rounds (default `200`).

- learning_rate:

  Shrinkage (default `0.1`).

- max_depth:

  Maximum tree depth (default `5`).

- pexp_bins:

  Number of piecewise-exponential time intervals (default `10`).

- ...:

  Additional arguments passed to
  [`fastgbm::fastgbm()`](https://rdrr.io/pkg/fastgbm/man/fastgbm.html)
  (for example `min_node_size`, `subsample`, `colsample`, `lambda`,
  `early_stopping` + `validation`, `seed`).

## Value

An object of class `"mlsurv_model"` with elements `model`, `learner`
(`"fastgbm"`), `formula`, `data`, `time`, `status`, and `objective`.

## See also

[`predict_fastgbm`](https://elbadisyimad.com/survalis/reference/predict_fastgbm.md),
[`fastgbm`](https://rdrr.io/pkg/fastgbm/man/fastgbm.html),
[`fit`](https://elbadisyimad.com/survalis/reference/fit.md) for the verb
interface

## Examples

``` r
mod <- fit_fastgbm(Surv(time, status) ~ age + karno + celltype, veteran,
                   ntrees = 50, pexp_bins = 8)
summary(mod)
#> 
#> ── fastgbm summary ─────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: fastgbm
#> Learner: fastgbm
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
