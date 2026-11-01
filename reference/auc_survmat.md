# Time-Dependent AUC from a Survival-Probability Matrix

Computes a cumulative/dynamic time-dependent AUC using predicted
survival probabilities at a specified time point (or the last column if
`t_star` is `NULL`), following the Uno et al. (2007) / Heagerty and
Zheng (2005) estimator as implemented in the timeROC package. Cases are
subjects with an observed event strictly before `t_star`; controls are
subjects known to survive beyond `t_star`. Subjects censored before
`t_star` are handled through IPCW weighting.

## Usage

``` r
auc_survmat(object, predicted, t_star = NULL)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object of
  length \\n\\.

- predicted:

  An `n x k` matrix or data frame of survival probabilities with columns
  named `"t=<time>"`.

- t_star:

  Optional numeric time at which to evaluate AUC; if omitted, the
  rightmost column of `predicted` is used.

## Value

A named numeric scalar: `"auc"`.

## Details

Risk scores are defined as `1 - S(t)` at the chosen time. The AUC is
computed over case-control pairs using inverse-probability-of-censoring
weights for cases and partial credit (0.5) for ties.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
sp <- matrix(
  stats::plogis(scale(veteran$karno)),
  ncol = 1,
  dimnames = list(NULL, "t=100")
)
auc_survmat(y, predicted = sp, t_star = 100)
#>      auc 
#> 0.795944 
```
