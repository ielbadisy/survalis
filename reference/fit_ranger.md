# Fit a Survival Random Forest Model Using ranger

This function fits a survival random forest model using the ranger
package and returns an object compatible with the `mlsurv_model` class.

## Usage

``` r
fit_ranger(formula, data, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the variables in the model.

- ...:

  Additional arguments passed to
  [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md).

## Value

An object of class `"mlsurv_model"` containing:

- `model` - the fitted `ranger` model

- `learner` - character string "ranger"

- `formula` - the model formula

- `data` - training data used to fit the model

## Details

This function wraps
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) for
survival analysis and stores the result in a standardized `mlsurv_model`
object.

## See also

[`predict_ranger`](https://elbadisyimad.com/survalis/reference/predict_ranger.md),
[`tune_ranger`](https://elbadisyimad.com/survalis/reference/tune_ranger.md),
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Examples

``` r
mod <- fit_ranger(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  num.trees = 25
)
summary(mod)
#> 
#> ── ranger summary ──────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: ranger
#> Learner: ranger
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
