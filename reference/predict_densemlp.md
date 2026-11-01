# Predict Survival Probabilities from a densemlp Model

Generates predicted survival probabilities at the requested `times` for
a model fitted with
[`fit_densemlp`](https://elbadisyimad.com/survalis/reference/fit_densemlp.md)
(`loss = "brier"`). The network returns survival on its own
discrete-time grid; values are linearly interpolated (and
constant-extrapolated) onto `times`.

## Usage

``` r
predict_densemlp(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` from
  [`fit_densemlp`](https://elbadisyimad.com/survalis/reference/fit_densemlp.md).

- newdata:

  A `data.frame` of predictors.

- times:

  Numeric vector of evaluation time points.

## Value

A `data.frame` of survival probabilities with one row per observation in
`newdata` and one column per time point (columns named `"t=<time>"`).

## See also

[`fit_densemlp`](https://elbadisyimad.com/survalis/reference/fit_densemlp.md)

## Examples

``` r
mod <- fit_densemlp(Surv(time, status) ~ age + karno, veteran,
                    n_bins = 8, epochs = 20)
predict_densemlp(mod, veteran[1:5, ], times = c(100, 200, 300))
#>        t=100      t=200      t=300
#> 1 0.07552381 0.03529123 0.01921588
#> 2 0.14958050 0.08408775 0.05788969
#> 3 0.10930761 0.07790300 0.07012705
#> 4 0.09201039 0.04241381 0.02772206
#> 5 0.14339453 0.07743385 0.05099604
```
