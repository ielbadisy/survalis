# Fit a Royston-Parmar Flexible Parametric Survival Model

Fits a Royston-Parmar flexible parametric survival model with the rpsurv
package and returns an object compatible with the `mlsurv_model`
interface. Restricted cubic splines model the baseline on the chosen
`scale`, optionally with time-varying covariate effects (`tve`), so
predicted survival curves are available at arbitrary time points without
a proportional-hazards assumption.

## Usage

``` r
fit_rpsurv(
  formula,
  data,
  df = 4,
  knots = NULL,
  tve = NULL,
  tve.df = 3,
  scale = c("hazard", "odds", "normal"),
  ...
)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the formula variables.

- df:

  Degrees of freedom for the baseline spline (default `4`).

- knots:

  Optional explicit knot vector (log time scale); default `NULL` derives
  knots from `df`.

- tve:

  Optional character vector of covariate names given time-varying
  effects via their own restricted cubic spline.

- tve.df:

  Degrees of freedom for each time-varying-effect spline (default `3`).

- scale:

  One of `"hazard"` (proportional hazards, default), `"odds"`, or
  `"normal"`.

- ...:

  Additional arguments passed to
  [`rpsurv::rpsurv()`](https://rdrr.io/pkg/rpsurv/man/rpsurv.html) (for
  example `control`).

## Value

An object of class `"mlsurv_model"` with elements `model`, `learner`
(`"rpsurv"`), `formula`, `data`, `time`, `status`, plus
`rhs_terms`/`cov_terms`/`xlev` used to rebuild the design matrix for
arbitrary `newdata` at predict time.

## See also

[`predict_rpsurv`](https://elbadisyimad.com/survalis/reference/predict_rpsurv.md),
[`rpsurv`](https://rdrr.io/pkg/rpsurv/man/rpsurv.html),
[`fit`](https://elbadisyimad.com/survalis/reference/fit.md) for the verb
interface

## Examples

``` r
mod <- fit_rpsurv(Surv(time, status) ~ age + karno + celltype, veteran, df = 3)
#> Warning: Hessian not invertible; standard errors unavailable
summary(mod)
#> 
#> ── rpsurv summary ──────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: rpsurv
#> Learner: rpsurv
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
