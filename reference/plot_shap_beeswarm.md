# SHAP Beeswarm Plot (Time-Integrated)

Plots the distribution of time-integrated SHAP contributions
([`compute_shap_matrix`](https://elbadisyimad.com/survalis/reference/compute_shap_matrix.md)
output) across observations, one row per feature, with points coloured
by each observation's (min-max scaled) raw feature value – summarizing
global feature impact across the entire prediction horizon in a single
panel.

## Usage

``` r
plot_shap_beeswarm(shap_matrix, top_n = NULL, title)
```

## Arguments

- shap_matrix:

  A data.frame returned by
  [`compute_shap_matrix`](https://elbadisyimad.com/survalis/reference/compute_shap_matrix.md).

- top_n:

  Optional integer; if supplied, only the `top_n` features with the
  largest mean absolute SHAP value are shown.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely.

## Value

A ggplot2 object.

## See also

[`compute_shap_matrix()`](https://elbadisyimad.com/survalis/reference/compute_shap_matrix.md)

## Examples

``` r
# \donttest{
mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
sm <- compute_shap_matrix(mod, newdata = veteran[1:15, ],
                          baseline_data = veteran, times = c(50, 100, 200),
                          sample.size = 5)
plot_shap_beeswarm(sm)

# }
```
