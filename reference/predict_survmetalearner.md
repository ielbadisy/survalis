# Predict with a Stacked Survival Meta‑Learner

Produces stacked survival probabilities by combining base learner
predictions via the time‑specific weights learned by
[`fit_survmetalearner`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md).

## Usage

``` r
predict_survmetalearner(model, newdata, times)
```

## Arguments

- model:

  A `"survmetalearner"` object returned by
  [`fit_survmetalearner()`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md).

- newdata:

  A data frame of new observations for prediction.

- times:

  Numeric vector of evaluation times. When they differ from the times
  the meta-learner was trained at, the weights are linearly interpolated
  in time.

## Value

A data frame with one row per observation and one column per requested
time (columns named `"t=<time>"`), containing stacked survival
probabilities.

## Details

For each base learner listed in `model$learners`, the corresponding
`predict_<learner>` function is called to obtain \\S\_\ell(t \mid x)\\.
Stacked predictions are computed as \\\hat{S}(t \mid x) = \sum\_\ell
w\_\ell(t) S\_\ell(t \mid x)\\, where \\w\_\ell(t)\\ are the learned
nonnegative weights for time \\t\\.

## See also

[`fit_survmetalearner`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md),
[`plot_survmetalearner_weights`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md)

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
predict_survmetalearner(meta_model, newdata = veteran[1:3, ], times = times)
#>        t=80     t=160
#> 1 0.4712827 0.1100416
#> 2 0.6155851 0.3162953
#> 3 0.4839684 0.2018965
```
