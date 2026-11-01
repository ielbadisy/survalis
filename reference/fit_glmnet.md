# Fit a Penalized Cox Proportional Hazards Model (glmnet)

Fits a Cox proportional hazards model with elastic net regularization
using glmnet. The function wraps
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html) to
select the optimal penalty parameter \\\lambda\\ by cross-validation.

## Usage

``` r
fit_glmnet(formula, data, alpha = 1, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- alpha:

  Numeric value in `[0, 1]` specifying the elastic net mixing parameter:

  - `1` for LASSO penalty

  - `0` for ridge penalty

  - between 0 and 1 for elastic net

- ...:

  Additional arguments passed to
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

## Value

An object of class `"mlsurv_model"` containing:

- `model` – the fitted `cv.glmnet` object

- `learner` – the string `"glmnet"`

- `formula`, `data`, `time`, and `status` metadata

## See also

[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
[`predict_glmnet`](https://elbadisyimad.com/survalis/reference/predict_glmnet.md)

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

summary(mod_glmnet$model)
#>            Length Class  Mode     
#> lambda     51     -none- numeric  
#> cvm        51     -none- numeric  
#> cvsd       51     -none- numeric  
#> cvup       51     -none- numeric  
#> cvlo       51     -none- numeric  
#> nzero      51     -none- numeric  
#> call        5     -none- call     
#> name        1     -none- character
#> glmnet.fit 12     coxnet list     
#> lambda.min  1     -none- numeric  
#> lambda.1se  1     -none- numeric  
#> index       2     -none- numeric  
```
