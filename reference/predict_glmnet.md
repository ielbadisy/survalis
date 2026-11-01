# Predict Survival Probabilities from a Penalized Cox Model (glmnet)

Predicts survival probabilities at specified time points from a fitted
penalized Cox proportional hazards model produced by
[`fit_glmnet`](https://elbadisyimad.com/survalis/reference/fit_glmnet.md).

## Usage

``` r
predict_glmnet(object, newdata, times, ...)
```

## Arguments

- object:

  A fitted `"mlsurv_model"` object from
  [`fit_glmnet`](https://elbadisyimad.com/survalis/reference/fit_glmnet.md).

- newdata:

  A `data.frame` of new observations for prediction.

- times:

  Numeric vector of time points at which to estimate survival
  probabilities. Must be in the same scale as the training time
  variable.

- ...:

  Not used.

## Value

A `data.frame` with one row per observation in `newdata` and one column
per requested time point (`"t=<time>"`).

## Details

Predictions are computed by:

1.  Computing the linear predictors for training and new data at
    `s = "lambda.min"`.

2.  Fitting a Cox PH model on the training linear predictors to estimate
    the baseline cumulative hazard.

3.  Interpolating the baseline hazard at each `times` and transforming
    via \\S(t \mid x) = \exp(-H_0(t) \exp(\eta))\\.

## See also

[`fit_glmnet`](https://elbadisyimad.com/survalis/reference/fit_glmnet.md)

## Examples

``` r

mod_glmnet <- fit_glmnet(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran
)
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

predict_glmnet(
  mod_glmnet,
  newdata = veteran[1:5, ],
  times = c(100, 200, 300)
)
#>       t=100     t=200      t=300
#> 1 0.4012275 0.1468623 0.06154467
#> 2 0.4984197 0.2316272 0.11934044
#> 3 0.4012275 0.1468623 0.06154467
#> 4 0.4012275 0.1468623 0.06154467
#> 5 0.4984197 0.2316272 0.11934044
```
