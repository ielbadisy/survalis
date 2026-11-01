# Convert a survival-probability matrix (survmat) to cumulative hazard

Uses the identity \\\Lambda(t \mid x) = -\log S(t \mid x)\\ to compute
cumulative hazards for each observation and time point.

## Usage

``` r
survmat_to_chf(S, eps = 1e-12)
```

## Arguments

- S:

  A `surv_mat`: numeric matrix/data.frame of survival probabilities.

- eps:

  Numeric in (0,1). Stabilizer to avoid `log(0)`. Values of `S` are
  clamped to \\\[\varepsilon, 1-\varepsilon\]\\ before logging. Default
  is `1e-12`.

## Value

A numeric matrix with the same dimensions as `S`, containing cumulative
hazards \\\Lambda(t \mid x)\\.

## Examples

``` r
S <- data.frame(`t=1` = c(0.9, 0.8), `t=2` = c(0.7, 0.6))
survmat_to_chf(S)
#>            t.1       t.2
#> [1,] 0.1053605 0.3566749
#> [2,] 0.2231436 0.5108256
```
