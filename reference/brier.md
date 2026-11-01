# Brier Score with IPCW for a Single Time Point

Computes the inverse-probability-of-censoring weighted (IPCW) Brier
score at a single time `t_star`.

## Usage

``` r
brier(object, pre_sp, t_star)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.

- pre_sp:

  Numeric vector of predicted survival probabilities \\S(t^{\*}\mid
  x_i)\\.

- t_star:

  Numeric evaluation time.

## Value

A named numeric scalar: `"brier"`.

## Details

The censoring distribution \\G(t)\\ is estimated via Kaplan-Meier on
`1 - status`. Observed events before `t_star` contribute \\S(t_i)^2 /
G(t_i)\\; those at risk at `t_star` contribute \\(1 - S(t^{\*}))^2 /
G(t^{\*})\\. Returns `NA` if \\G(t^{\*})\\ is undefined or zero.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
pre_sp <- stats::plogis(scale(veteran$karno))
brier(y, pre_sp = pre_sp, t_star = 100)
#>    brier 
#> 0.189903 
```
