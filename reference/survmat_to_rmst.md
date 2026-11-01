# Compute restricted mean survival time (RMST) from a survival-probability matrix (survmat)

Computes \\\mathrm{RMST}(\tau) = \int_0^\tau S(t)\\ dt\\ for each
observation, via the trapezoidal rule on the provided time grid
(delegated to
[`tvrmst::rmst_dynamic()`](https://rdrr.io/pkg/tvrmst/man/rmst_dynamic.html)),
followed by linear interpolation to `tau` when `tau` falls between grid
points.

## Usage

``` r
survmat_to_rmst(S, times, tau = max(times))
```

## Arguments

- S:

  A `surv_mat`: numeric matrix/data.frame of survival probabilities.

- times:

  Numeric vector of time points corresponding to columns of `S`. Must be
  strictly increasing and non-negative.

- tau:

  Positive numeric scalar; upper integration limit. Default is
  `max(times)`.

## Value

Numeric vector of RMST values (length = `nrow(S)`).

## Details

If the grid does not include `0`, the function prepends \\(0, S(0)=1)\\
to correctly integrate from time 0. When `tau` exceeds the last grid
point, the value is held at the definite integral over the full grid (no
extrapolation past the observed curve).

## Examples

``` r
times <- c(1, 2, 3)
S <- matrix(c(0.9, 0.8, 0.7,
              0.95,0.9,0.85), nrow = 2, byrow = TRUE)
colnames(S) <- paste0("t=", times)
survmat_to_rmst(S, times, tau = 3)
#> [1] 2.550 2.775
```
