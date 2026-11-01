# Predict Survival Probabilities from an `rpart` Survival Tree

Predicts survival probabilities at specified time points from a fitted
`rpart` survival tree model, assuming an exponential survival
distribution within each terminal node.

## Usage

``` r
predict_rpart(object, newdata, times, ...)
```

## Arguments

- object:

  An `"mlsurv_model"` object created by
  [`fit_rpart()`](https://elbadisyimad.com/survalis/reference/fit_rpart.md).

- newdata:

  A `data.frame` containing the predictor variables for prediction.

- times:

  Numeric vector of time points at which to compute survival
  probabilities.

- ...:

  Additional arguments passed to internal prediction functions.

## Value

A `data.frame` with one row per observation in `newdata` and one column
per requested time point, containing predicted survival probabilities.

## Details

Predictions are based on converting the predicted mean survival times
from the survival tree into survival probabilities under an exponential
assumption: \$\$S(t) = \exp(-t / \hat{\mu})\$\$ where \\\hat{\mu}\\ is
the predicted mean survival time from the terminal node.

## Examples

``` r

mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200      t=300
#> 1 0.4062035 0.1650013 0.06702412
#> 2 0.6080413 0.3697142 0.22480154
#> 3 0.4062035 0.1650013 0.06702412
#> 4 0.4062035 0.1650013 0.06702412
#> 5 0.6080413 0.3697142 0.22480154
```
