# Score a Fitted Survival Model on Its Training Data

Computes one or more performance metrics for a fitted `mlsurv_model`,
predicting on the training data with the model's corresponding
`predict_*` function.

## Usage

``` r
score_survmodel(
  model,
  times,
  metrics = c("cindex", "ibs", "brier", "iae", "ise")
)
```

## Arguments

- model:

  An object of class `"mlsurv_model"`.

- times:

  Numeric vector of evaluation times. For `"brier"`, must be a single
  time.

- metrics:

  Character vector of metrics to compute. Supported: `"cindex"`,
  `"auc"`, `"ibs"`, `"brier"`, `"iae"`, `"ise"`, `"ece"`.

## Value

A data.frame with columns `metric` and `value`.

## Details

The function constructs the appropriate `predict_*` function name from
`model$learner`, predicts survival probabilities on `model$data`, builds
a `Surv` object from `model$formula`, and computes the metrics. If
`"brier"` is requested with multiple `times`, an error is thrown.

## Examples

``` r
fitted_model <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
score_survmodel(
  fitted_model,
  times = c(80, 160),
  metrics = c("cindex", "ibs")
)
#>        metric value
#> cindex cindex 0.711
#> ibs       ibs 0.176
```
