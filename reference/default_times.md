# Default evaluation-time grid

Builds the time grid used by the verb-level functions
([`predict()`](https://rdrr.io/r/stats/predict.html),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md))
when the caller does not pass `times` explicitly. The grid spans the
observed event times so that integrated metrics (IBS, IAE, ISE) and
survival-curve summaries are computed over the region where the data are
informative.

## Usage

``` r
default_times(time, status = NULL, n = 50L, range = c(0.05, 0.95))
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- status:

  Optional 0/1 event indicator aligned with `time`. When supplied, the
  grid is built from event times only; otherwise all positive times are
  used.

- n:

  Integer; maximum number of grid points (default `50`).

- range:

  Two-value probability vector giving the lower and upper quantiles of
  the event-time distribution spanned by the grid (default
  `c(0.05, 0.95)`), which trims unstable tails.

## Value

A numeric vector, strictly increasing and positive, of length at most
`n`.

## Details

The chosen grid is recorded on every result object (`res$times`) and
printed, so that runs remain reproducible and comparable even though the
default is data-dependent.

## Examples

``` r
default_times(veteran$time, veteran$status)
#>  [1]   7.000000   8.000000   8.030612  10.347959  12.000000  13.026531
#>  [7]  15.345918  18.000000  19.000000  20.000000  21.000000  24.000000
#> [13]  25.000000  26.348980  30.000000  31.000000  34.344898  42.005102
#> [19]  44.337755  48.670408  51.000000  52.000000  52.668367  54.002041
#> [25]  59.667347  68.996939  79.992857  82.663265  88.992857  94.990816
#> [31]  99.000000 101.986735 109.974490 111.327551 117.000000 121.971429
#> [37] 132.325510 139.658163 143.990816 153.970408 162.000000 176.854082
#> [43] 200.321429 223.848980 241.854082 265.748980 285.608163 339.602041
#> [49] 379.904082 390.950000
default_times(veteran$time, veteran$status, n = 10)
#>  [1]   7.00  15.00  23.50  35.45  52.15  83.70 111.00 145.75 230.85 390.95
```
