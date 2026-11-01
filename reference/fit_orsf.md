# Fit an Oblique Random Survival Forest (ORSF) Model

Fits an Oblique Random Survival Forest using the aorsf package. This
method builds an ensemble of oblique decision trees for survival
analysis, where splits are based on linear combinations of features,
allowing for improved performance in high-dimensional or correlated
feature settings.

## Usage

``` r
fit_orsf(formula, data, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A data frame containing the variables specified in `formula`.

- ...:

  Additional arguments passed to
  [`aorsf::orsf()`](https://docs.ropensci.org/aorsf/reference/orsf.html).

## Value

An object of class `"mlsurv_model"` containing:

- `model`: The fitted aorsf ORSF model object.

- `learner`: Character string, always `"orsf"`.

- `formula`: The survival formula used.

- `data`: The training dataset.

- `time`: Name of the survival time variable.

- `status`: Name of the event indicator variable.

## Details

ORSF models extend traditional Random Survival Forests by allowing
oblique splits, which can improve prediction accuracy when predictors
are correlated. Missing data are omitted by default
(`na_action = "omit"`).

## References

Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI,
Bandyopadhyay D. Oblique random survival forests. *Annals of Applied
Statistics*. 2019;13(3):1847-1883. doi:10.1214/19-AOAS1261

## Examples

``` r
mod <- fit_orsf(Surv(time, status) ~ age + karno, data = veteran)
summary(mod)
#> 
#> ── orsf summary ────────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno
#> Engine: aorsf
#> Learner: orsf
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
