# Predict Survival Probabilities from a Survival SVM

Turns the fitted log-time predictor \\f(x)\\ of a model from
[`fit_survsvm()`](https://elbadisyimad.com/survalis/reference/fit_survsvm.md)
into survival probabilities at the requested times.

## Usage

``` r
predict_survsvm(
  object,
  newdata,
  times,
  dist = c("km", "lognormal", "exp", "weibull"),
  shape = 1
)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_survsvm()`](https://elbadisyimad.com/survalis/reference/fit_survsvm.md).

- newdata:

  A `data.frame` of new observations with the predictors used in the
  formula.

- times:

  Numeric vector of evaluation time points (same scale as the training
  times).

- dist:

  How the distribution of the error \\\log T - f(x)\\ is obtained:

  `"km"`

  :   (default) Kaplan-Meier estimate of the cross-fitted residuals,
      \\S(t \mid x) = \hat S\_\varepsilon(\log t - f(x))\\.

  `"lognormal"`

  :   normal errors, mean and scale estimated by censored maximum
      likelihood on the cross-fitted residuals.

  `"exp"`, `"weibull"`

  :   the earlier fixed mappings with \\e^{f(x)}\\ as the time scale:
      \\\exp(-t / e^{f(x)})\\ and \\\exp(-(t / e^{f(x)})^{shape})\\.

- shape:

  Weibull shape parameter (only for `dist = "weibull"`).

## Value

A base `data.frame` with one row per observation in `newdata` and
columns `t=<time>` holding survival probabilities, non-increasing in
time.

## See also

[`fit_survsvm()`](https://elbadisyimad.com/survalis/reference/fit_survsvm.md),
[`tune_survsvm()`](https://elbadisyimad.com/survalis/reference/tune_survsvm.md)

## Examples

``` r
mod <- fit_survsvm(Surv(time, status) ~ age + celltype + karno, data = veteran)
predict_survsvm(mod, newdata = veteran[1:5, ], times = c(100, 300, 500))
#>       t=100      t=300      t=500
#> 1 0.4911426 0.11209330 0.04670554
#> 2 0.5367487 0.14887392 0.06538776
#> 3 0.4450980 0.09341109 0.04670554
#> 4 0.4988167 0.11209330 0.04670554
#> 5 0.5367487 0.14887392 0.06538776
predict_survsvm(mod, newdata = veteran[1:5, ], times = c(100, 300, 500),
                dist = "lognormal")
#>       t=100     t=300      t=500
#> 1 0.4309795 0.1401000 0.06666084
#> 2 0.4887962 0.1751322 0.08765499
#> 3 0.3982848 0.1222569 0.05647869
#> 4 0.4363685 0.1431729 0.06845023
#> 5 0.4862965 0.1735202 0.08666106
```
