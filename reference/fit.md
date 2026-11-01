# Fit a survival learner

The entry point of the survalis verb interface. `fit()` selects a
learner by id, passes engine arguments through `spec`, and returns a
`survalis_fit` object that the other verbs
([`predict()`](https://rdrr.io/r/stats/predict.html),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md),
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md))
consume.

## Usage

``` r
fit(formula, data, model, spec = list(), seed = NULL, ...)
```

## Arguments

- formula:

  A survival formula with a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) outcome on the
  left-hand side, e.g. `Surv(time, status) ~ age + karno`.

- data:

  A data frame containing the formula variables.

- model:

  A learner id from `list_survlearners()$learner` (for example
  `"coxph"`, `"ranger"`, `"densemlp"`).

- spec:

  A named list of engine arguments forwarded to the underlying
  `fit_<model>()` (for example `list(num.trees = 500)` for `"ranger"`).

- seed:

  Optional integer seed set immediately before fitting, for learners
  with a stochastic training step.

- ...:

  Not used; engine arguments passed here give an error. Use `spec`.

## Value

An object of class `c("survalis_fit", "mlsurv_model")` with elements
`model`, `learner`, `engine`, `formula`, `data`, `time`, `status`,
`model_id`, and `spec`.

## Details

`survalis_fit` extends the internal `mlsurv_model` contract, so every
legacy `predict_*()` / `compute_*()` helper keeps working on it
unchanged.

## See also

[`predict.survalis_fit()`](https://elbadisyimad.com/survalis/reference/predict.survalis_fit.md),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md),
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md),
[`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md)

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno + celltype, veteran,
         model = "ranger", spec = list(num.trees = 200))
m
#> <survalis_fit> model = ranger  (engine: ranger)
#>   formula: Surv(time, status) ~ age + karno + celltype
#>   data:    137 obs, 128 events (93%)
#>   spec:    num.trees = 200
```
