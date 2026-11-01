# The survmat prediction contract

Internal constructor and validator for `survmat`, the single prediction
container used by the verb-level
[`predict()`](https://rdrr.io/r/stats/predict.html) method. A `survmat`
is an `n x length(times)` numeric matrix carrying its evaluation grid in
the `"times"` attribute (column names `"t=<time>"` are kept for display
and for backward compatibility with the `*_survmat` metric helpers) and
the predicted quantity in the `"quantity"` attribute.

## Usage

``` r
new_survmat(x, times, quantity = "survival")

validate_survmat(x)

as_survmat(x, times, quantity = "survival")
```

## Arguments

- x:

  A numeric matrix (or data frame / vector coercible to one) of
  predictions, one row per observation and one column per entry of
  `times`.

- times:

  Numeric vector, strictly increasing and non-negative, of length
  `ncol(x)`.

- quantity:

  One of `"survival"`, `"risk"`, `"chf"`, `"hazard"`. Controls printing
  and the `[0, 1]` range check (applied to `"survival"` and `"risk"`
  only).

## Value

A `survmat`: a numeric matrix with class
`c("survmat", "matrix", "array")`, attribute `"times"`, and attribute
`"quantity"`.

## Details

`new_survmat()` builds the object without checks; `validate_survmat()`
enforces the contract; `as_survmat()` is the checked constructor used
internally after a learner's raw `predict_*()` output has been
finalised.
