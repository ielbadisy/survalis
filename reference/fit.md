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
  `"coxph"`, `"ranger"`, `"densemlp"`), or a custom learner: a list with
  `fit = function(formula, data, ...)` returning any fitted object and
  `predict = function(object, newdata, times)` returning an
  `nrow(newdata) x length(times)` matrix of survival probabilities. An
  optional `name` element labels it in printed output. A custom learner
  works with every verb.

- spec:

  A named list of engine arguments forwarded to the learner's fitting
  routine (for example `list(num.trees = 500)` for `"ranger"`).

- seed:

  Optional integer seed set immediately before fitting, for learners
  with a stochastic training step.

- ...:

  Not used; engine arguments passed here give an error. Use `spec`.

## Value

An object of class `c("survalis_fit", "mlsurv_model")` with elements
`model`, `learner`, `engine`, `formula`, `data`, `time`, `status`,
`model_id`, and `spec`.

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

# A custom learner: an exponential AFT model from survival::survreg
expo <- list(
  name = "exponential",
  fit = function(formula, data, ...) survival::survreg(formula, data, dist = "exponential"),
  predict = function(object, newdata, times) {
    rate <- 1 / predict(object, newdata, type = "response")
    exp(-outer(rate, times))
  }
)
m2 <- fit(Surv(time, status) ~ age + karno, veteran, model = expo)
predict(m2, veteran[1:3, ], times = c(90, 180))
#> <survmat> survival  |  3 obs x 2 times  [90, 180]
#>           t=90     t=180
#> [1,] 0.4587666 0.2104668
#> [2,] 0.5774304 0.3334259
#> [3,] 0.4625431 0.2139461
```
