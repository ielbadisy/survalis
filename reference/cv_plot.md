# Boxplot of Cross-Validation Metric Distributions

Visualizes per-fold metric values from
[`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
using a boxplot with jittered points.

## Usage

``` r
cv_plot(cv_results, title)
```

## Arguments

- cv_results:

  A tibble/data frame as returned by
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md).

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

## Examples

``` r
cv_results <- data.frame(
  metric = c("cindex", "cindex", "ibs", "ibs"),
  value = c(0.62, 0.66, 0.19, 0.21)
)
cv_plot(cv_results)
```
