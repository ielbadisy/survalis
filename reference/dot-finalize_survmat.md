# Standardize learner prediction output to the survmat contract

Internal utility used by `predict_*()` methods and evaluators to coerce
raw predictions to a validated survival-probability matrix contract.

## Usage

``` r
.finalize_survmat(S, times, clamp = TRUE, enforce_monotone = TRUE)
```

## Arguments

- S:

  Raw survival predictions as matrix, data.frame, or vector.

- times:

  Numeric vector of evaluation times aligned with columns of `S`.

- clamp:

  Logical; if `TRUE`, clamp values to `[0, 1]`.

- enforce_monotone:

  Logical; if `TRUE`, enforce non-increasing survival over increasing
  `times` via cumulative minima.

## Value

A base data.frame with columns named `t=<time>`.
