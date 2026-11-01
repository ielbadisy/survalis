# Predict Survival Probabilities from an ORSF Model

Generates survival probability predictions at specified time points from
a fitted Oblique Random Survival Forest model.

## Usage

``` r
predict_orsf(object, newdata, times, ...)
```

## Arguments

- object:

  A fitted ORSF model object from
  [`fit_orsf`](https://elbadisyimad.com/survalis/reference/fit_orsf.md).

- newdata:

  A data frame of new observations for prediction.

- times:

  A numeric vector of time points at which to estimate survival
  probabilities.

- ...:

  Additional arguments passed to `predict.aorsf()`.

## Value

A data frame where each row corresponds to an observation in `newdata`
and each column corresponds to a requested prediction time
(`"t=<time>"`).

## Details

Predictions are computed using the `pred_type = "surv"` option from
aorsf, returning estimated survival probabilities for each observation
at the specified time points.

## Examples

``` r
mod <- fit_orsf(Surv(time, status) ~ age + karno, data = veteran)
pred <- predict_orsf(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
head(pred)
#>       t=100     t=200      t=300
#> 1 0.2961354 0.1374661 0.05922410
#> 2 0.5919963 0.2334583 0.13539055
#> 3 0.5525881 0.3405800 0.17234644
#> 4 0.5158287 0.2265305 0.07801944
#> 5 0.5899023 0.2232810 0.11557720
```
