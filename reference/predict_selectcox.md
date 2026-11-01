# Predict Survival Probabilities with a Selected Cox Model

Generates survival probabilities at specified time points using a fitted
selected‑predictor Cox model (`mlsurv_model`) via
[`pec::predictSurvProb()`](https://rdrr.io/pkg/pec/man/predictSurvProb.html).

## Usage

``` r
predict_selectcox(object, newdata, times)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_selectcox()`](https://elbadisyimad.com/survalis/reference/fit_selectcox.md).

- newdata:

  A `data.frame` of new observations for prediction. Must include the
  same predictor variables used at training time.

- times:

  Numeric vector of evaluation time points (same scale as the training
  survival time). Must be non‑negative and finite.

## Value

A base `data.frame` with one row per observation in `newdata` and
columns named `t={time}` (character), containing numeric survival
probabilities in \[0, 1\].

## Details

Internally calls `pec::predictSurvProb(object$model, newdata, times)`
and renames columns to the standard `t=...` convention used by
`survalis`.

## See also

[`fit_selectcox()`](https://elbadisyimad.com/survalis/reference/fit_selectcox.md),
[`tune_selectcox()`](https://elbadisyimad.com/survalis/reference/tune_selectcox.md)

## Examples

``` r
mod <- fit_selectcox(Surv(time, status) ~ age + celltype + karno, data = veteran)
predict_selectcox(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.5933344 0.3296536 0.2021624
#> 2 0.6820581 0.4433239 0.3097821
#> 3 0.5933344 0.3296536 0.2021624
#> 4 0.5933344 0.3296536 0.2021624
#> 5 0.6820581 0.4433239 0.3097821
```
