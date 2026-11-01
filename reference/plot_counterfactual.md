# Plot Counterfactual Recommendations

Visualizes counterfactual feature changes returned by
[`compute_counterfactual`](https://elbadisyimad.com/survalis/reference/compute_counterfactual.md),
ranking recommendations by either raw survival gain or penalized gain.

## Usage

``` r
plot_counterfactual(
  counterfactual_df,
  metric = c("penalized_gain", "survival_gain"),
  top_n = NULL,
  include_negative = FALSE,
  title
)
```

## Arguments

- counterfactual_df:

  A data frame returned by
  [`compute_counterfactual()`](https://elbadisyimad.com/survalis/reference/compute_counterfactual.md).

- metric:

  Character scalar; one of `"penalized_gain"` (default) or
  `"survival_gain"` indicating which column to plot on the x-axis.

- top_n:

  Optional integer limiting the plot to the top `n` features after
  sorting by `metric`.

- include_negative:

  Logical; if `FALSE` (default), only positive-gain recommendations are
  shown.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

## Examples

``` r
df <- veteran
df$A <- df$trt
mod <- fit_coxph(survival::Surv(time, status) ~ A + age + karno, data = df)
cf <- compute_counterfactual(
  model = mod,
  newdata = df[1, , drop = FALSE],
  times = c(50, 100, 150),
  target_time = 100,
  features_to_change = c("A", "age", "karno"),
  grid.size = 10
)
plot_counterfactual(cf)
```
