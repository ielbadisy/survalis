# Validate and standardize a survival-probability matrix (survmat)

Internal utility to validate the basic structure of a `surv_mat` and
coerce it to a numeric matrix. Optionally checks column names against
`times`.

## Usage

``` r
.survmat_as_matrix(S, times = NULL, strict_colnames = FALSE, clamp = FALSE)
```

## Arguments

- S:

  A `surv_mat`: numeric matrix or data.frame of survival probabilities.

- times:

  Optional numeric vector of time points. If provided, requires
  `length(times) == ncol(S)` and enforces the column names `"t=<time>"`.

- strict_colnames:

  Logical; if `TRUE` and `times` is provided, checks that existing
  colnames are exactly `"t=<time>"` (after coercion to character). If
  `FALSE`, colnames are overwritten to match the contract.

- clamp:

  Logical; if `TRUE`, clamps `S` into `[0,1]`.

## Value

A numeric matrix.
