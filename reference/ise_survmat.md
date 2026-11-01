# Integrated Squared Error Against Kaplan-Meier

Computes \\\int ( \bar{S}(t) - \hat{S}\_{KM}(t) )^2 \\ dt\\ where
\\\bar{S}(t)\\ is the mean predicted survival and \\\hat{S}\_{KM}(t)\\
is the Kaplan-Meier curve. Integration uses KM event times and a left
Riemann sum.

## Usage

``` r
ise_survmat(object, sp_matrix, times)
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

A named numeric scalar: `"ise"`.

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
times <- c(60, 120)
lp <- stats::plogis(scale(veteran$karno))
sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
colnames(sp) <- c("t=60", "t=120")
ise_survmat(y, sp_matrix = sp, times = times)
#>      ise 
#> 139.7427 
```
