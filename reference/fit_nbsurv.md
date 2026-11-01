# Fit a Conditional Naive Bayes Survival Model

Fits a conditional naive Bayes survival model with the nbsurv package
and returns an object compatible with the `mlsurv_model` interface.
Frames survival at each horizon as a binary classification problem
("survive beyond t" vs. "fail before t") and applies Bayes' rule under a
conditional-independence assumption over predictors, with IPCW weighting
for censoring.

## Usage

``` r
fit_nbsurv(formula, data, scale = TRUE, laplace = 1, min_sd = 0.05, ...)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the formula variables.

- scale:

  Logical; standardize continuous predictors using training mean/sd
  (default `TRUE`).

- laplace:

  Laplace smoothing constant for categorical predictors (default `1`).

- min_sd:

  Floor for estimated SDs of Gaussian class-conditional distributions
  (default `0.05`).

- ...:

  Additional arguments passed to
  [`nbsurv::nbsurv()`](https://rdrr.io/pkg/nbsurv/man/nbsurv.html) (for
  example `time_grid`, `cov_structure`, `shrinkage`, `time_smooth`,
  `bandwidth`).

## Value

An object of class `"mlsurv_model"` with elements `model`, `learner`
(`"nbsurv"`), `formula`, `data`, `time`, and `status`.

## See also

[`predict_nbsurv`](https://elbadisyimad.com/survalis/reference/predict_nbsurv.md),
[`nbsurv`](https://rdrr.io/pkg/nbsurv/man/nbsurv.html),
[`fit`](https://elbadisyimad.com/survalis/reference/fit.md) for the verb
interface

## Examples

``` r
mod <- fit_nbsurv(Surv(time, status) ~ age + karno + celltype, veteran)
summary(mod)
#> 
#> ── nbsurv summary ──────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: nbsurv
#> Learner: nbsurv
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
