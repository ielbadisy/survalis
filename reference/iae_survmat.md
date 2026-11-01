# Integrated Absolute Error Against Kaplan-Meier

Computes \\\int \mid \bar{S}(t) - \hat{S}\_{KM}(t) \mid \\ dt\\ where
\\\bar{S}(t)\\ is the mean predicted survival across subjects and
\\\hat{S}\_{KM}(t)\\ is the Kaplan-Meier estimate. Integration is
carried out over KM event times using a left Riemann sum.

## Usage

``` r
iae_survmat(object, sp_matrix, times)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.

- sp_matrix:

  Matrix/data frame of survival probabilities (rows = subjects, columns
  aligned with `times`).

- times:

  Numeric vector of times corresponding to the columns of `sp_matrix`.

## Value

A named numeric scalar: `"iae"`.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
times <- c(60, 120)
lp <- stats::plogis(scale(veteran$karno))
sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
colnames(sp) <- c("t=60", "t=120")
iae_survmat(y, sp_matrix = sp, times = times)
#>      iae 
#> 355.5755 
```
