# Fit an Accelerated Failure Time Model Using Generalized Estimating Equations

Fits an accelerated failure time (AFT) model via aftgee, which
implements GEE methodology for censored survival data. Returns an
`mlsurv_model`-compatible object used throughout the package.

## Usage

``` r
fit_aftgee(formula, data, corstr = "independence")
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A data.frame containing the variables in the model.

- corstr:

  Working correlation structure; one of `"independence"`,
  `"exchangeable"`, or `"ar1"`. Default: `"independence"`.

## Value

An object of class `mlsurv_model` with components:

- `model`: the fitted `aftgee` model.

- `learner`: `"aftgee"`.

- `formula`: the survival formula.

- `data`: the training data.

## Details

The AFT model assumes \\\log(T) = X \beta + \epsilon\\, where \\T\\ is
survival time, \\X\\ is the covariate matrix, and \\\epsilon\\ is an
error term whose distribution determines the baseline survival. aftgee
estimates \\\beta\\ using GEE. In this wrapper we set `id = NULL`,
assuming one observation per subject.

## References

Jin, Z., Lin, D. Y., Wei, L. J., & Ying, Z. (2003). Rank-based inference
for the accelerated failure time model. *Biometrika*, 90(2), 341-353.

## Examples

``` r
# \donttest{
mod <- fit_aftgee(
  Surv(time, status) ~ trt + karno + age,
  data = veteran
)
head(predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120)))
#>        t=20      t=60     t=120
#> 1 0.9207975 0.6224206 0.3514893
#> 2 0.9608786 0.7461304 0.4877205
#> 3 0.8557272 0.4851251 0.2325601
#> 4 0.9103448 0.5964888 0.3267598
#> 5 0.9618224 0.7497250 0.4922122
# }
```
