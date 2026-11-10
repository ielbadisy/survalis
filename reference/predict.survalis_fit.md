# Predict from a fitted survival learner

Produces predictions from a
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) object.
Every matrix-valued `type` returns a `survmat` (an `n x length(times)`
numeric matrix carrying its grid in the `"times"` attribute); the
scalar-valued types (`"rmst"`, `"quantile"`, `"median"`) return a
numeric vector (or an `n x length(probs)` matrix for a vector `probs`).

## Usage

``` r
# S3 method for class 'survalis_fit'
predict(
  object,
  newdata = NULL,
  times = NULL,
  type = c("survival", "risk", "chf", "hazard", "rmst", "quantile", "median"),
  ...
)
```

## Arguments

- object:

  A `survalis_fit` from
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md).

- newdata:

  A data frame of predictors. Defaults to the training data. Rows with
  missing values in model variables are rejected.

- times:

  Numeric evaluation grid. When `NULL`, up to 50 event-time quantiles
  between the 5th and 95th percentiles are used and recorded on the
  result.

- type:

  One of:

  `"survival"`

  :   predicted survival probability \\S(t)\\ (default).

  `"risk"`

  :   cumulative risk \\1 - S(t)\\.

  `"chf"`

  :   cumulative hazard \\-\log S(t)\\.

  `"hazard"`

  :   piecewise-constant hazard on the grid.

  `"rmst"`

  :   restricted mean survival time; upper limit `tau` (`...`), default
      `max(times)`.

  `"quantile"`

  :   survival-time quantile(s); requires `probs` (`...`).

  `"median"`

  :   median survival time (`probs = 0.5`).

- ...:

  Type-specific arguments: `tau` for `"rmst"`, `probs` for `"quantile"`.

## Value

A `survmat` for `"survival"`, `"risk"`, `"chf"`, `"hazard"`; a numeric
vector for `"rmst"`, `"median"` and scalar-`probs` `"quantile"`; an
`n x length(probs)` matrix for vector `probs`.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno, veteran, "ranger",
         spec = list(num.trees = 100))
predict(m, veteran[1:5, ], times = c(90, 180, 365))
#> <survmat> survival  |  5 obs x 3 times  [90, 365]
#>           t=90     t=180     t=365
#> [1,] 0.2624219 0.1967211 0.1828112
#> [2,] 0.8763340 0.4313011 0.3130361
#> [3,] 0.7686504 0.4792258 0.2617366
#> [4,] 0.7931211 0.3105448 0.1516631
#> [5,] 0.8492272 0.3460189 0.2651582
predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "risk")
predict(m, veteran[1:5, ], type = "rmst", tau = 365)
#>        1        2        3        4        5 
#> 119.4606 212.7497 181.8384 171.1826 188.9202 
```
