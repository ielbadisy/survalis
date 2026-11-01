# Resolve a caller-supplied `times` argument against the default grid

Resolve a caller-supplied `times` argument against the default grid

## Usage

``` r
.resolve_times(times, time, status = NULL, n = 50L)
```

## Arguments

- times:

  `NULL` or a numeric vector.

- time, status:

  Passed to `default_times()` when `times` is `NULL`.

- n:

  Grid size for the default.

## Value

A validated numeric vector (strictly increasing, positive).
