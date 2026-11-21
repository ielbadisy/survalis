# Plot an evaluation

`type = "metrics"` draws one panel per metric with the per-split values
as a box plot and points; each panel has its own value axis and its
strip says whether higher or lower is better. `type = "time"` draws the
time-specific metrics (`"brier"` evaluated at several `times`, stored as
`"brier@<time>"`) against time: the mean over the splits with its 95%
interval, or the single value when the predictions were scored once.

## Usage

``` r
# S3 method for class 'survalis_eval'
plot(
  x,
  type = c("metrics", "time"),
  reference = NULL,
  labels = NULL,
  xlab = "Time",
  time.labels = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ...
)
```

## Arguments

- x:

  A `survalis_eval` object.

- type:

  `"metrics"` (default) or `"time"`.

- reference:

  For `type = "time"`: an optional second `survalis_eval` on the same
  times, drawn dashed in grey, for example the Kaplan-Meier curve given
  to everyone
  ([`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  on a matrix).

- labels:

  For `type = "time"`: legend labels of `x` and `reference`. `NULL` uses
  their model names.

- xlab:

  For `type = "time"`: the time axis label.

- time.labels:

  For `type = "time"`: `NULL` labels the evaluation times as they are; a
  function maps them to labels, e.g. `function(t) t / 365`.

- title, subtitle, caption:

  Plot title, subtitle and caption; `NULL` (default) draws none.

- ...:

  Unused.

## Value

A ggplot.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
ev <- evaluate(m, times = c(30, 90, 180), resampling = cv(3, seed = 1),
               metrics = c("cindex", "brier"))
plot(ev)

km <- summary(survival::survfit(Surv(time, status) ~ 1, veteran), times = c(30, 90, 180))$surv
ref <- evaluate(matrix(km, nrow(veteran), 3, byrow = TRUE),
                y = Surv(veteran$time, veteran$status), times = c(30, 90, 180),
                metrics = "brier")
plot(ev, type = "time", reference = ref, labels = c("Cox model", "Kaplan-Meier"),
     xlab = "Days")
```
