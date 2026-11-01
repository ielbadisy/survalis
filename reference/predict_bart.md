# Predict Survival Probabilities from a BART Survival Model

Generates survival probability predictions at requested times for new
data using a model fitted by
[`fit_bart`](https://elbadisyimad.com/survalis/reference/fit_bart.md).

## Usage

``` r
predict_bart(object, newdata, times)
```

## Arguments

- object:

  An `mlsurv_model` returned by
  [`fit_bart`](https://elbadisyimad.com/survalis/reference/fit_bart.md).

- newdata:

  A `data.frame` of new observations for prediction.

- times:

  Numeric vector of time points at which to estimate survival
  probabilities (same scale as the training time).

## Value

A base `data.frame` with one row per observation in `newdata` and
columns named `"t=<time>"` (character), containing survival
probabilities in \[0, 1\].

## Details

The BART engine predicts survival on its internal grid
`object$eval_times`. Requested `times` are aligned to that grid by
nearest‐neighbor matching, returning one survival estimate per requested
time.

## See also

[`fit_bart`](https://elbadisyimad.com/survalis/reference/fit_bart.md),
[`surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html)

## Examples

``` r
ex_data <- veteran[1:40, c("time", "status", "age", "karno", "celltype")]
mod_bart <- fit_bart(
  Surv(time, status) ~ age + karno + celltype,
  data = ex_data,
  K = 1,
  ntree = 5,
  ndpost = 20,
  nskip = 5,
  mc.cores = 1,
  seed = 42
)
predict_bart(mod_bart, newdata = ex_data[1:5, ], times = c(10, 30, 60))
#>           t=10         t=30         t=60
#> 1 0.1600936718 0.1600936718 0.1600936718
#> 2 0.0268880343 0.0268880343 0.0268880343
#> 3 0.0039681229 0.0039681229 0.0039681229
#> 4 0.0009865883 0.0009865883 0.0009865883
#> 5 0.0003031816 0.0003031816 0.0003031816
```
