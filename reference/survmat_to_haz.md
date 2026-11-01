# Convert a survival-probability matrix (survmat) to hazards on a time grid

Computes *individual-specific* (i.e., conditional on covariates) hazards
from predicted survival curves \\S(t \mid x_i)\\ on a discrete time
grid.

## Usage

``` r
survmat_to_haz(S, times, eps = 1e-12, t0 = 0)
```

## Arguments

- S:

  A `surv_mat`: numeric matrix/data.frame of survival probabilities.

- times:

  Numeric vector of time points corresponding to the columns of `S`.
  Must be strictly increasing and non-negative.

- eps:

  Numeric in (0,1). Stabilizer to avoid `log(0)` and division by zero;
  survival probabilities are clamped to \\\[\varepsilon,
  1-\varepsilon\]\\. Default is `1e-12`.

- t0:

  Numeric scalar; left boundary for the first interval. Default is `0`.
  If your grid already starts at 0 and you want the first interval to be
  (0, t1\], keep `t0 = 0`.

## Value

A numeric matrix of hazards with the same dimensions as `S`. Each column
corresponds to the interval ending at `times[j]`.

## Details

This function returns a *piecewise-constant hazard* per interval: \$\$
h\_{i,j} \approx \frac{\Lambda_i(t_j) - \Lambda_i(t\_{j-1})}{t_j -
t\_{j-1}} = \frac{\log S_i(t\_{j-1}) - \log S_i(t_j)}{t_j - t\_{j-1}}.
\$\$ with \\S_i(0) = 1\\ and \\\Lambda_i(t) = -\log S_i(t)\\.

## Examples

``` r
times <- c(1, 2, 3)
S <- matrix(c(0.9, 0.8, 0.7,
              0.95,0.9,0.85), nrow = 2, byrow = TRUE)
colnames(S) <- paste0("t=", times)
H <- survmat_to_haz(S, times)
H
#>                        t=1        t=2
#> [1,] 0.10536052 0.11778304 0.13353139
#> [2,] 0.05129329 0.05406722 0.05715841
```
