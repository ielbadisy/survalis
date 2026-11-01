# Predict Survival from an Aalen Additive Hazards Model

Computes survival probabilities at specified time points from a model
fitted with
[`fit_aalen`](https://elbadisyimad.com/survalis/reference/fit_aalen.md).

## Usage

``` r
predict_aalen(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` returned by
  [`fit_aalen`](https://elbadisyimad.com/survalis/reference/fit_aalen.md).

- newdata:

  A data frame of new observations.

- times:

  Numeric vector of time points at which to evaluate survival
  probabilities.

## Value

A data frame (rows = observations, columns = `"t=<time>"`).

## Examples

``` r
# \donttest{
  mod <- fit_aalen(Surv(time, status) ~ trt + karno + age, data = veteran, max.time = 600)
  head(predict_aalen(mod, newdata = veteran[1:5, ], times = 0:10))
#>   t=0       t=1       t=2       t=3       t=4       t=5       t=6      t=7
#> 1   1 0.9974340 0.9974340 0.9974340 0.9974340 0.9974340 0.9974340 0.997434
#> 2   1 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.000000
#> 3   1 1.0000000 0.9531478 0.9190050 0.8750720 0.8750720 0.8750720 0.875072
#> 4   1 0.9990563 0.9990563 0.9990563 0.9990563 0.9990563 0.9990563 0.982176
#> 5   1 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.000000
#>         t=8       t=9      t=10
#> 1 0.9895864 0.9608233 0.9547647
#> 2 1.0000000 0.9919142 0.9904614
#> 3 0.8720192 0.8720192 0.8502718
#> 4 0.9656562 0.9446462 0.9335840
#> 5 1.0000000 0.9947253 0.9941717
# }
```
