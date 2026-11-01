# Predict Survival Probabilities from a flexsurvreg Model

This function generates survival probability predictions at specified
time points from a model fitted with
[`fit_flexsurvreg`](https://elbadisyimad.com/survalis/reference/fit_flexsurvreg.md).

## Usage

``` r
predict_flexsurvreg(object, newdata, times, ...)
```

## Arguments

- object:

  A fitted `mlsurv_model` object returned by
  [`fit_flexsurvreg`](https://elbadisyimad.com/survalis/reference/fit_flexsurvreg.md).

- newdata:

  A data frame of new observations for which to predict survival
  probabilities.

- times:

  Numeric vector of time points at which to estimate survival
  probabilities.

- ...:

  Additional arguments passed to
  [`predict.flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/predict.flexsurvreg.md).

## Value

A data frame of survival probabilities with one row per observation in
`newdata` and one column per time point. Column names are of the form
`"t=100"`, `"t=200"`, etc.

## Examples

``` r
mod_flex <- fit_flexsurvreg(Surv(time, status) ~ age + celltype + karno,
                            data = veteran,
                            dist = "weibull")
predict_flexsurvreg(mod_flex, newdata = veteran[1:5, ],
                    times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.6148791 0.3610455 0.2080287
#> 2 0.6964750 0.4687296 0.3110514
#> 3 0.5721321 0.3104626 0.1648523
#> 4 0.6068279 0.3512133 0.1993623
#> 5 0.6975962 0.4703117 0.3126710
```
