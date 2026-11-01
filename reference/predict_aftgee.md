# Predict Survival Probabilities from an `aftgee` Model

Computes survival probabilities \\S(t \mid x)\\ at specified time points
from a fitted AFT model using a log-normal approximation for the error
distribution.

## Usage

``` r
predict_aftgee(object, newdata, times = NULL)
```

## Arguments

- object:

  An `mlsurv_model` produced by
  [`fit_aftgee`](https://elbadisyimad.com/survalis/reference/fit_aftgee.md).

- newdata:

  A data.frame of predictor values for prediction.

- times:

  Numeric vector of time points at which to estimate survival
  probabilities. If `NULL`, a sequence of 20 values from `1` to the
  maximum observed time in `object$data` is used.

## Value

A data.frame of survival probabilities with one row per `newdata`
observation and one column per requested time (`"t=<time>"`).

## Details

We use the approximation \$\$S(t \mid x) \approx 1 -
\Phi\\\left(\frac{\log t - x^\top \hat{\beta}}{\sigma}\right),\$\$ where
\\\Phi\\ is the standard normal CDF. Here we set \\\sigma = 1\\ as a
simple, distribution-agnostic proxy; this yields a monotone-in-time
score useful for benchmarking. For production use, prefer a parametric
AFT fit where \\\sigma\\ is estimated.

## Examples

``` r
  mod <- fit_aftgee(Surv(time, status) ~ trt + karno + age, data = veteran)
  predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120))
#>        t=20      t=60     t=120
#> 1 0.9207975 0.6224206 0.3514893
#> 2 0.9608786 0.7461304 0.4877205
#> 3 0.8557272 0.4851251 0.2325601
#> 4 0.9103448 0.5964888 0.3267598
#> 5 0.9618224 0.7497250 0.4922122
```
