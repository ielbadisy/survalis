# Decision-weighted scores for survival predictions

Scores predicted risks at one or more horizons with a decision-weighted
proper score: the loss of the rule "act if predicted risk exceeds `c`",
averaged over thresholds `c` with a weight `w(c)`. The weight states
which risk thresholds matter. Every such score equals a weighted area
under the net benefit curve, with weight `(1 - c) w(c)`, so a threshold
range restricts model comparison to the decisions clinicians actually
make.

## Usage

``` r
dwscore(
  object,
  y,
  times = NULL,
  weight = c("range", "survacc", "brier", "log", "beta"),
  thresholds = c(0.05, 0.2),
  shape = c(2, 18)
)
```

## Arguments

- object:

  A numeric matrix of survival probabilities (rows = subjects, columns =
  times), for example from
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- y:

  A [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object with one entry per row of `object`.

- times:

  Evaluation times. `NULL` reads the `"times"` attribute or column names
  of the form `"t=<time>"`.

- weight:

  Threshold weight:

  `"range"`

  :   uniform on `thresholds` (default); only decisions with a threshold
      in that range count.

  `"survacc"`

  :   spherical weight \\(c^2 + (1 - c)^2)^{-3/2}\\; the loss is
      `1 - SurvAcc`.

  `"brier"`

  :   constant weight 2; the loss is the Brier score.

  `"log"`

  :   weight \\1 / (c (1 - c))\\; the loss is the logarithmic score.

  `"beta"`

  :   Beta density with parameters `shape`.

- thresholds:

  Length-2 numeric range `c(lower, upper)` in `[0, 1]` for
  `weight = "range"` (default `c(0.05, 0.20)`).

- shape:

  Length-2 numeric Beta parameters for `weight = "beta"` (default
  `c(2, 18)`, mean threshold 0.10).

## Value

A data frame with columns `time`, `weight` and `loss` (lower is better
for every weight).

## Details

With risk \\r = 1 - S(t \mid X)\\ and event status \\D = I(T \le t)\\,
the loss is \\L_w(r, 0) = \int_0^r w(c)\\ c\\ dc\\ for a subject event
free at `t` and \\L_w(r, 1) = \int_r^1 w(c) (1 - c)\\ dc\\ for a subject
with an event by `t`. Right censoring is handled by inverse probability
of censoring weights from the Kaplan-Meier estimate of the censoring
distribution: subjects event free at `t` are weighted by \\1 / G(t)\\,
subjects with an observed event at \\U \le t\\ by \\1 / G(U-)\\, and
subjects censored before `t` get weight zero. This assumes censoring is
independent of the covariates. Range and Beta weights integrate to one,
so those losses lie on the scale of an average regret per decision.

In
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
and [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md), use
the metric ids `"survacc"` (SurvAcc, higher is better) and
`"range:<lower>-<upper>"`, for example `"range:0.05-0.20"` (loss, lower
is better). Both need a single time.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno, veteran[1:100, ], "coxph")
S <- predict(m, veteran[-(1:100), ], times = c(90, 180))
y <- Surv(veteran$time[-(1:100)], veteran$status[-(1:100)])
dwscore(S, y)
#>   time          weight       loss
#> 1   90 range(0.05,0.2) 0.03659910
#> 2  180 range(0.05,0.2) 0.01463964
dwscore(S, y, weight = "brier")
#>   time weight       loss
#> 1   90  brier 0.17687881
#> 2  180  brier 0.09568638
dwscore(S, y, weight = "range", thresholds = c(0.10, 0.30))
#>   time         weight       loss
#> 1   90 range(0.1,0.3) 0.05470029
#> 2  180 range(0.1,0.3) 0.02342342

# the same scores inside the verbs (single time)
S180 <- predict(m, veteran[-(1:100), ], times = 180)
evaluate(S180, y = y, metrics = c("brier", "survacc", "range:0.05-0.20"))
#> <survalis_eval> model = supplied predictions
#>   scored once on 37 subjects (no resampling)
#>   grid: 1 times in [180, 180]
#>   brier           0.0960
#>   survacc         0.8980
#>   range:0.05-0.20 0.0150
```
