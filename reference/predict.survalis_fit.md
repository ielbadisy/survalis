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
  Arguments of the learner's own prediction function are passed on, e.g.
  `dist` and `shape` for `"survsvm"`.

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
#> [1,] 0.2910391 0.2251554 0.2047503
#> [2,] 0.9093085 0.4370865 0.3389879
#> [3,] 0.8037139 0.4674085 0.2577191
#> [4,] 0.7375774 0.3164041 0.1609650
#> [5,] 0.8691375 0.3258824 0.2514397
predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "risk")
#> <survmat> risk  |  5 obs x 3 times  [90, 365]
#>            t=90     t=180     t=365
#> [1,] 0.70896089 0.7748446 0.7952497
#> [2,] 0.09069146 0.5629135 0.6610121
#> [3,] 0.19628609 0.5325915 0.7422809
#> [4,] 0.26242265 0.6835959 0.8390350
#> [5,] 0.13086253 0.6741176 0.7485603
predict(m, veteran[1:5, ], type = "rmst", tau = 365)
#>        1        2        3        4        5 
#> 126.9180 218.6376 182.2984 168.1643 185.9864 
```
