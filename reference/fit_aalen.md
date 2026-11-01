# Fit an Additive Hazards (Aalen) Model

Fits an additive hazards regression model using timereg's
[`aalen`](https://rdrr.io/pkg/timereg/man/aalen.html) and returns an
`mlsurv_model` compatible with the survalis workflow.

## Usage

``` r
fit_aalen(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A data frame containing the variables in `formula`.

- max.time:

  Optional maximum follow-up time used by the fitting routine.

- n.sim:

  Integer; number of simulations for variance estimation (default `0`).

- resample.iid:

  Integer; indicator for iid resampling (passed to `aalen`).

## Value

A list of class `"mlsurv_model"` with elements: `model`,
`learner="aalen"`, `engine="timereg"`, `formula`, `data`, `time`,
`status`.

## Details

The Aalen model assumes an additive hazard: \$\$\lambda(t \mid X) =
\beta_0(t) + X^\top \beta(t),\$\$ with nonparametric cumulative
coefficient functions.

## Examples

``` r
# \donttest{
  mod_aalen <- fit_aalen(
    Surv(time, status) ~ trt + karno + age,
    data = veteran
  )
  head(predict_aalen(mod_aalen, newdata = veteran[1:5, ], times = c(50, 100, 150)))
#>        t=50     t=100     t=150
#> 1 0.6663730 0.4568039 0.2844424
#> 2 0.8555688 0.6748876 0.3966400
#> 3 0.5254661 0.3474935 0.2589745
#> 4 0.6364268 0.4332509 0.2793249
#> 5 0.8621505 0.6808684 0.3978419
# }
```
