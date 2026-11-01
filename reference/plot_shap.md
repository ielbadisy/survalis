# Plot SHAP-like contributions for survival models

Plots time-dependent SHAP estimates (lines over time) or aggregated SHAP
(bar chart) returned by
[`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md).

## Usage

``` r
plot_shap(shapley_result, type = c("auto"), title)
```

## Arguments

- shapley_result:

  A data frame returned by
  [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md).
  If it has a `time` column, a time plot is produced; otherwise an
  aggregated bar plot.

- type:

  One of `"auto"` (default), `"time"`, or `"integrated"`. With `"auto"`,
  the presence of a `time` column determines the plot type.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
shap_td <- compute_shap(mod, veteran[10, , drop = FALSE],
veteran, times = c(50, 100), sample.size = 5)
p1 <- plot_shap(shap_td)      # auto -> time plot

shap_ag <- compute_shap(mod, veteran[10, , drop = FALSE],
veteran, times = c(50, 100),
sample.size = 5, aggregate = TRUE, method = "meanabs")
p2 <- plot_shap(shap_ag)      # auto -> bar plot
# }
```
