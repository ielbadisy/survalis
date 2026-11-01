# Evaluation grid of a survmat

Returns the numeric time grid of a `survmat`, reading the `"times"`
attribute and falling back to parsing `"t=<time>"` column names for
objects produced by the legacy `predict_*()` path.

## Usage

``` r
survmat_times(x)
```

## Arguments

- x:

  A `survmat` or a legacy survival-probability matrix / data frame.

## Value

A numeric vector.
