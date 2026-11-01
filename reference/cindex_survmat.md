# Concordance Index from a Survival-Probability Matrix

Computes Harrell's concordance index using predicted survival
probabilities at a specified time point (or the last column if `t_star`
is `NULL`).

## Usage

``` r
cindex_survmat(object, predicted, t_star = NULL)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object of
  length \\n\\.

- predicted:

  An `n x k` matrix or data frame of survival probabilities with columns
  named `"t=<time>"`.

- t_star:

  Optional numeric time at which to evaluate the c-index; if omitted,
  the rightmost column of `predicted` is used.

## Value

A named numeric scalar: `"C index"`.

## Details

Risk scores are defined as `1 - S(t)` at the chosen time. Ties receive
partial credit (0.5). Pairs not comparable due to censoring are
excluded.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
sp <- matrix(
  stats::plogis(scale(veteran$karno)),
  ncol = 1,
  dimnames = list(NULL, "t=100")
)
cindex_survmat(y, predicted = sp, t_star = 100)
#> C index 
#> 0.70847 
```
