# Integrated Brier Score (Discrete Integration)

Computes the Integrated Brier Score (IBS) over a vector of times, using
discrete left Riemann integration of IPCW Brier scores.

## Usage

``` r
ibs_survmat(object, sp_matrix, times)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.

- sp_matrix:

  Matrix/data frame of survival probabilities with one column per time
  in `times`.

- times:

  Numeric vector of strictly increasing times (length must equal
  `ncol(sp_matrix)`).

## Value

A named numeric scalar: `"ibs"`.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
times <- c(60, 120)
lp <- stats::plogis(scale(veteran$karno))
sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
colnames(sp) <- c("t=60", "t=120")
ibs_survmat(y, sp_matrix = sp, times = times)
#>      ibs 
#> 0.174546 
```
