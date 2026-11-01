# Fit a Componentwise Gradient Boosted Cox Model (blackboost)

Fits a survival boosting model using mboost's `blackboost` with a Cox
proportional hazards loss
([`mboost::CoxPH()`](https://rdrr.io/pkg/mboost/man/Family.html)) and
shallow tree base-learners (partykit). Returns an `mlsurv_model`
compatible with the `survalis` pipeline.

## Usage

``` r
fit_blackboost(
  formula,
  data,
  weights = NULL,
  mstop = 100,
  nu = 0.1,
  minsplit = 10,
  minbucket = 4,
  maxdepth = 2,
  ...
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- weights:

  Optional case weights passed to
  [`mboost::blackboost()`](https://rdrr.io/pkg/mboost/man/blackboost.html).

- mstop:

  Integer; number of boosting iterations (stopping iteration). Default
  `100`.

- nu:

  Learning rate/shrinkage in `(0,1]`; default `0.1`.

- minsplit, :

  minbucket, maxdepth Tree control parameters passed via
  [`partykit::ctree_control()`](https://rdrr.io/pkg/partykit/man/ctree_control.html)
  (defaults: `10`, `4`, `2`).

- minbucket:

  Minimum bucket size per tree node.

- maxdepth:

  Maximum tree depth.

- ...:

  Additional arguments forwarded to
  [`mboost::boost_control()`](https://rdrr.io/pkg/mboost/man/control.html)
  and/or
  [`mboost::blackboost()`](https://rdrr.io/pkg/mboost/man/blackboost.html).

## Value

An object of class `mlsurv_model` with elements:

- model:

  Fitted `mboost` model.

- learner:

  `"blackboost"`.

- formula, data, time, status:

  Original inputs/metadata.

## Details

The base-learner is a conditional inference tree controlled by
[`partykit::ctree_control()`](https://rdrr.io/pkg/partykit/man/ctree_control.html).
The loss is the partial likelihood for Cox PH via
[`mboost::CoxPH()`](https://rdrr.io/pkg/mboost/man/Family.html). Use
`mstop` and `nu` to control complexity.

## References

Bühlmann P, Hothorn T (2007). Boosting algorithms: regularization,
prediction and model fitting. *Statistical Science*.

## Examples

``` r
# \donttest{
  mod <- fit_blackboost(
    Surv(time, status) ~ age + karno + celltype,
    data = veteran
  )
  head(predict_blackboost(mod, newdata = veteran[1:5, ], times = c(100, 200)))
#>       t=100     t=200
#> 1 0.4849286 0.2374386
#> 2 0.5369659 0.2907356
#> 3 0.4849286 0.2374386
#> 4 0.4849286 0.2374386
#> 5 0.5369659 0.2907356
# }
```
