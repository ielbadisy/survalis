# Fit a Random Survival Forest (RSF) Model

Fits a Random Survival Forest for right-censored time-to-event data
using randomForestSRC and returns a standardized `mlsurv_model`
compatible with the `survalis` framework.

## Usage

``` r
fit_rsf(formula, data, ntree = 500, mtry = NULL, nodesize = 15, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the variables referenced in `formula`.

- ntree:

  Integer; number of trees to grow (default: 500).

- mtry:

  Integer or `NULL`; number of variables randomly selected at each
  split. If `NULL`, the engine's default is used.

- nodesize:

  Integer; minimum terminal node size (default: 15).

- ...:

  Additional arguments forwarded to
  [`rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Value

An object of class `mlsurv_model`, a named list with elements:

- model:

  The fitted
  [`randomForestSRC::rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html)
  object.

- learner:

  Character scalar identifying the learner (`"rsf"`).

- engine:

  Character scalar naming the engine (`"randomForestSRC"`).

- formula:

  The original survival formula.

- data:

  The training dataset (or a minimal subset needed for prediction).

- time:

  Name of the survival time variable.

- status:

  Name of the event indicator (1 = event, 0 = censored).

## Details

RSF extends random forests to survival data by growing an ensemble of
survival trees on bootstrap samples and aggregating survival functions
across trees.

## References

Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS (2008). Random survival
forests. *Annals of Applied Statistics*, 2(3):841-860.

## Examples

``` r
mod_rsf <- fit_rsf(Surv(time, status) ~ age + celltype + karno,
                   data = veteran, ntree = 200)

times <- c(100, 200, 300)
pred_probs <- predict_rsf(mod_rsf, newdata = veteran[1:5, ], times = times)
print(round(pred_probs, 3))
#>   t=100 t=200 t=300
#> 1 0.621 0.319 0.138
#> 2 0.824 0.496 0.298
#> 3 0.663 0.408 0.210
#> 4 0.730 0.453 0.244
#> 5 0.811 0.464 0.272
```
