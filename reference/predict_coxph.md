# Predict Survival Probabilities from a Cox PH Model

Generates predicted survival probabilities at specified time points for
new data, using a fitted Cox proportional hazards model.

## Usage

``` r
predict_coxph(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` object returned by
  [`fit_coxph()`](https://elbadisyimad.com/survalis/reference/fit_coxph.md).

- newdata:

  A data frame containing the predictor variables for which to compute
  predictions.

- times:

  A numeric vector of time points at which to evaluate survival
  probabilities.

## Value

A data frame of survival probabilities with columns named `"t=<time>"`
for each requested time.

## Details

Predictions are computed using
[`pec::predictSurvProb()`](https://rdrr.io/pkg/pec/man/predictSurvProb.html).
The output is formatted as a data frame with one row per observation in
`newdata` and one column per time point.

## Examples

``` r
mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
predict_coxph(mod_cox, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.6142681 0.3541697 0.2175837
#> 2 0.6944383 0.4599242 0.3194203
#> 3 0.5556797 0.2860796 0.1589937
#> 4 0.6033305 0.3408724 0.2056867
#> 5 0.6959633 0.4620783 0.3216209
```
