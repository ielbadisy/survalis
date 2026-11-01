# Tidy data frames from verb result objects

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
that return the most useful rectangular view of each result object, for
downstream plotting or reporting:

## Usage

``` r
# S3 method for class 'survalis_eval'
as.data.frame(x, ...)

# S3 method for class 'survalis_tune'
as.data.frame(x, ...)

# S3 method for class 'survalis_compare'
as.data.frame(x, ...)

# S3 method for class 'survalis_estimate'
as.data.frame(x, ...)
```

## Arguments

- x:

  A verb result object.

- ...:

  Ignored.

## Value

A data frame.

## Details

- `survalis_eval` – the per-split metric values (`id`, `metric`,
  `value`).

- `survalis_tune` – the scored hyperparameter grid.

- `survalis_compare` – the model-by-metric leaderboard (`model`,
  `metric`, `mean`, `se`).

- `survalis_estimate` – a one-row summary (`estimand`, `tau`,
  `estimate`, `conf.low`, `conf.high`, `model`).
