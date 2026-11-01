# Plot Time‑Varying Stacking Weights

Visualizes the learned nonnegative NNLS stacking weights \\w\_\ell(t)\\
over time for each base learner.

## Usage

``` r
plot_survmetalearner_weights(model, title)
```

## Arguments

- model:

  A `"survmetalearner"` object from
  [`fit_survmetalearner`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md).

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object showing weight trajectories (one line per learner).

## Examples

``` r
form <- Surv(time, status) ~ age + karno + trt
times <- c(80, 160)
mod_cox <- fit_coxph(form, data = veteran)
mod_rpart <- fit_rpart(form, data = veteran)
base_models <- list(coxph = mod_cox, rpart = mod_rpart)
base_preds <- list(
  coxph = predict_coxph(mod_cox, newdata = veteran, times = times),
  rpart = predict_rpart(mod_rpart, newdata = veteran, times = times)
)
meta_model <- fit_survmetalearner(
  base_preds = base_preds,
  time = veteran$time,
  status = veteran$status,
  times = times,
  base_models = base_models,
  formula = form,
  data = veteran
)
plot_survmetalearner_weights(meta_model)
```
