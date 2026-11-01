# Infer time points from survmat column names

Internal helper that parses `"t=<time>"` column names when an explicit
time vector is not supplied.

## Usage

``` r
.infer_survmat_times(S)
```

## Arguments

- S:

  A survival-probability matrix or data frame.

## Value

Numeric vector of inferred times.
