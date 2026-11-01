# Time-Dependent ROC/AUC Curve from a Survival-Probability Matrix

Computes the cumulative/dynamic time-dependent AUC (Uno et al., 2007;
Heagerty and Zheng, 2005) at each time in `times`, using marginal
(Kaplan-Meier) inverse-probability-of-censoring weighting for cases.
This is a thin vectorized wrapper around
[`auc_survmat`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)
that returns a full \\AUC(t)\\ curve instead of a single time point,
matching the quantity reported by
[`timeROC::timeROC()`](https://rdrr.io/pkg/timeROC/man/timeROC.html)
with `weighting = "marginal"`.

## Usage

``` r
timeroc_survmat(object, predicted, times)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object of
  length \\n\\.

- predicted:

  An `n x k` matrix or data frame of survival probabilities with columns
  named `"t=<time>"`, one column per entry of `times`.

- times:

  Numeric vector of evaluation time points; each must have a matching
  `"t=<time>"` column in `predicted`.

## Value

A `data.frame` (class `"timeroc_survmat"`) with columns `time` and
`auc`, one row per entry of `times`.

## References

Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for
t-year survivors with censored regression models.
[doi:10.1198/016214507000000149](https://doi.org/10.1198/016214507000000149)

Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
curves.
[doi:10.1111/j.0006-341X.2005.030814.x](https://doi.org/10.1111/j.0006-341X.2005.030814.x)

## See also

[`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
times <- c(30, 90, 180)
sp <- sapply(times, function(t) stats::plogis(-scale(veteran$karno) + log(t) / 5))
colnames(sp) <- paste0("t=", times)
timeroc_survmat(y, predicted = sp, times = times)
#>   time      auc
#> 1   30 0.142451
#> 2   90 0.172116
#> 3  180 0.287582
```
