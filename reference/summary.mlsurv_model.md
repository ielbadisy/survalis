# Summarize an `mlsurv_model`

Produces a compact, human‑readable overview of a fitted survival learner
that follows the `mlsurv_model` contract (e.g., objects returned by
`fit_*()` in this toolkit). The summary prints the learner id, engine,
original formula, and basic data characteristics (sample size, predictor
names, time range, event rate). A structured list is returned invisibly
for programmatic use.

## Usage

``` r
# S3 method for class 'mlsurv_model'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"mlsurv_model"` produced by a `fit_*()` function.

- ...:

  Ignored; included for S3 signature compatibility.

## Value

Invisibly returns a list of class `"summary.mlsurv_model"` containing:

- learner:

  Character id of the learner (e.g., `"ranger"`, `"coxph"`).

- engine:

  Underlying package/engine used to fit the model.

- formula:

  The original survival formula.

- data_summary:

  List with `observations`, `predictors`, `time_range`, and
  `event_rate`.

## Details

This method relies on the presence of the fields `learner`, `formula`,
and `data` stored in the fitted object (the standard `mlsurv_model`
contract). Output printing uses cli if available.

## See also

[`fit_coxph()`](https://elbadisyimad.com/survalis/reference/fit_coxph.md),
other `fit_*()` learners returning `mlsurv_model` objects.

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + trt + celltype, data = veteran)
summary(mod)
#> 
#> ── coxph summary ───────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + trt + celltype
#> Engine: survival
#> Learner: coxph
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, trt, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"

s <- summary(mod)  # capture the structured result invisibly
#> 
#> ── coxph summary ───────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + trt + celltype
#> Engine: survival
#> Learner: coxph
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, trt, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
str(s)
#> List of 4
#>  $ learner     : chr "coxph"
#>  $ engine      : chr "survival"
#>  $ formula     :Class 'formula'  language Surv(time, status) ~ age + trt + celltype
#>   .. ..- attr(*, ".Environment")=<environment: 0x636e3071c0c0> 
#>  $ data_summary:List of 4
#>   ..$ observations: int 137
#>   ..$ predictors  : chr [1:5] "age" "trt" "celltypesmallcell" "celltypeadeno" ...
#>   ..$ time_range  : num [1:2] 1 999
#>   ..$ event_rate  : num 0.934
#>  - attr(*, "class")= chr "summary.mlsurv_model"
```
