# Fit a Cox Proportional Hazards Model

Fits a Cox proportional hazards regression model using the
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
function, and returns an object compatible with the `mlsurv_model`
interface.

## Usage

``` r
fit_coxph(formula, data, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A data frame containing the variables in the model.

- ...:

  Additional arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

An object of class `"mlsurv_model"` containing:

- `model` – the fitted `coxph` object

- `learner` – the string `"coxph"`

- `formula` – the survival formula used

- `data` – the training dataset

- `time`, `status` – names of the survival outcome variables

## Details

The fitted object is stored along with metadata such as the learner name
(`"coxph"`), original formula, data, and names of the time and status
variables. The function requires the survival and pec packages.

## Examples

``` r
mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
summary(mod_cox)
#> 
#> ── coxph summary ───────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: survival
#> Learner: coxph
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
