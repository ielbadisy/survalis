# Plot Interaction Strengths for Survival Models

Visualizes interaction outputs from
[`compute_interactions`](https://elbadisyimad.com/survalis/reference/compute_interactions.md)
as (i) a ranked bar chart for one-way interactions, (ii) a pairwise
heatmap, or (iii) time-varying interaction trajectories.

## Usage

``` r
plot_interactions(object, type = c("1way", "heatmap", "time"), title)
```

## Arguments

- object:

  A data frame returned by
  [`compute_interactions()`](https://elbadisyimad.com/survalis/reference/compute_interactions.md)
  whose columns match the requested `type`.

- type:

  One of `"1way"`, `"heatmap"`, or `"time"`.

- title:

  Plot title. If missing, an automatically generated title is used (none
  for `"1way"`). Pass `NULL` to omit the title entirely (e.g., for
  journals requiring caption-only figures).

## Value

A ggplot2 object.

## Details

**1way**: Bars rank features by Friedman-H interaction strength at the
target time. **heatmap**: Tiles show pairwise interaction magnitudes
(symmetric). **time**: Lines show interaction strength vs. time for each
feature.

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
times <- c(80, 160)
ia <- compute_interactions(
  model = mod,
  data = veteran,
  times = times,
  target_time = 80,
  features = c("age", "karno"),
  type = "1way",
  grid.size = 6
)
plot_interactions(ia, type = "1way")

```
