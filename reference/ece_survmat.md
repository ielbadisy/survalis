# Expected Calibration Error (ECE) at a Single Time Point

Computes a binned Expected Calibration Error at time `t_star`, using
quantile-based bins of predicted survival probabilities.

## Usage

``` r
ece_survmat(object, sp_matrix, t_star, n_bins = 10L, p = 1, weighted = TRUE)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.

- sp_matrix:

  Matrix or data frame of survival probabilities with rows = subjects
  and (optionally) columns named `"t=<time>"`.

- t_star:

  Numeric evaluation time (must be a single value).

- n_bins:

  Integer number of quantile bins (default `10`).

- p:

  Power for the L^p error (default `1` gives absolute error).

- weighted:

  Logical; if `TRUE` (default) weights bin errors by `n_k / N` where
  `n_k` is the bin size.

## Value

A named numeric scalar: `"ece"`.

## Details

Let \\\hat{S}\_k(t^\*)\\ be the mean predicted survival in bin \\k\\ and
\\\tilde{S}\_k(t^\*)\\ the Kaplan-Meier survival at \\t^\*\\ in that
bin. The (weighted) ECE is \$\$ \mathrm{ECE}(t^\*) = \left( \sum_k w_k
\\ \|\hat{S}\_k(t^\*) - \tilde{S}\_k(t^\*)\|^p \right)^{1/p} \$\$ with
\\w_k = n_k / N\\.

## See also

The calibration *curve* is available as
`interpret(fit, method = "calibration")` (via
[`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md));
`ece_survmat()` is the scalar summary of the same idea used by
`evaluate(metrics = "ece")`.

## Examples

``` r
y <- survival::Surv(
  time = c(1, 2, 3, 4, 6, 7, 8, 9),
  event = c(1, 1, 0, 1, 0, 1, 1, 0)
)
sp <- cbind("t=5" = c(0.15, 0.20, 0.35, 0.40, 0.55, 0.65, 0.75, 0.80))
ece_survmat(y, sp_matrix = sp, t_star = 5, n_bins = 4)
#>     ece 
#> 0.29375 
```
