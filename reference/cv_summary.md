# Summarize Cross-Validation Results

Produces mean, standard deviation, standard error, and 95\\ for each
metric returned by
[`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md).

## Usage

``` r
cv_summary(cv_results, digits = 3)
```

## Arguments

- cv_results:

  A tibble/data frame as returned by
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md).

- digits:

  Integer number of decimal places for `mean`, `sd`, `se`, `lower`,
  `upper` (default `3`).

## Value

A data.frame with columns: `metric`, `mean`, `sd`, `n`, `se`, `lower`,
`upper`.

## Examples

``` r
cv_results <- data.frame(
  metric = c("cindex", "cindex", "ibs", "ibs"),
  value = c(0.62, 0.66, 0.19, 0.21)
)
cv_summary(cv_results)
#>   metric mean    sd n   se lower upper
#> 1 cindex 0.64 0.028 2 0.02 0.601 0.679
#> 2    ibs 0.20 0.014 2 0.01 0.180 0.220
```
