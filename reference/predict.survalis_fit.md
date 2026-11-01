# Predict from a fitted survival learner

Produces predictions from a
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) object.
Every matrix-valued `type` returns a
[survmat](https://elbadisyimad.com/survalis/reference/survmat.md) (an
`n x length(times)` numeric matrix carrying its grid in the `"times"`
attribute); the scalar-valued types (`"rmst"`, `"quantile"`, `"median"`)
return a numeric vector (or an `n x length(probs)` matrix for a vector
`probs`).

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

  Numeric evaluation grid. When `NULL`,
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md)
  is used and recorded on the result.

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

A [survmat](https://elbadisyimad.com/survalis/reference/survmat.md) for
`"survival"`, `"risk"`, `"chf"`, `"hazard"`; a numeric vector for
`"rmst"`, `"median"` and scalar-`probs` `"quantile"`; an
`n x length(probs)` matrix for vector `probs`.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno, veteran, "ranger",
         spec = list(num.trees = 100))
predict(m, veteran[1:5, ], times = c(90, 180, 365))
#> <survmat> survival  |  5 obs x 3 times  [90, 365]
#>           t=90     t=180     t=365
#> [1,] 0.2710445 0.2053296 0.1857899
#> [2,] 0.9244360 0.4663988 0.3378462
#> [3,] 0.7278194 0.4639897 0.2474468
#> [4,] 0.7714770 0.3181532 0.1623141
#> [5,] 0.8235893 0.3335018 0.2486018
predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "risk")
predict(m, veteran[1:5, ], type = "rmst", tau = 365)
#>        1        2        3        4        5 
#> 120.0569 223.4384 172.5615 169.1469 183.1916 
```
