# Parse a survival outcome from a model formula

Internal helper that extracts the time and status components from a
formula whose left-hand side is `Surv(...)` or `pkg::Surv(...)`.

## Usage

``` r
.parse_surv_formula(formula, data)
```

## Arguments

- formula:

  A model formula with a survival outcome on the left-hand side.

- data:

  A data.frame used to resolve status recoding expressions such as
  `status == 2`.

## Value

A list with `time_col`, `status_col`, `event_value`, and
`recode_status`.
