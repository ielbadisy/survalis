# Compute a survival-time quantile from a survival-probability matrix (survmat) with grid-based approach

Returns the smallest grid time \\t\\ such that \\S(t) \le 1-p\\. This is
a right-continuous approximation on the provided time grid.

## Usage

``` r
survmat_to_quantile(S, times, p = 0.5)
```

## Arguments

- S:

  A `surv_mat`: numeric matrix/data.frame of survival probabilities.

- times:

  Numeric vector of time points corresponding to columns of `S`. Must be
  strictly increasing and non-negative.

- p:

  Probability in (0,1). Default is `0.5`.

## Value

Numeric vector of quantile times (length = `nrow(S)`). Returns `NA` if
the threshold is not crossed on the grid (e.g., survival stays above
\\1-p\\ over all `times`).

## Details

For `p = 0.5`, this yields a grid-based **median survival time**.

## Examples

``` r
times <- c(1, 2, 3)
S <- matrix(c(0.9, 0.6, 0.4,
              0.95,0.9,0.85), nrow = 2, byrow = TRUE)
colnames(S) <- paste0("t=", times)
survmat_to_quantile(S, times, p = 0.5)  # median on the grid
#> [1]  3 NA
```
