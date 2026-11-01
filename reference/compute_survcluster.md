# Cluster Predicted Survival Curves

Clusters the individualized survival curves predicted by a fitted
`mlsurv_model` using
[`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html)
(Partitioning Around Medoids on a weighted feature representation of the
curves).

## Usage

``` r
compute_survcluster(
  model,
  newdata,
  times,
  K = NULL,
  K_max = 10,
  distance = c("L2", "L1"),
  seed = NULL,
  ...
)
```

## Arguments

- model:

  An `mlsurv_model` fitted via a `fit_*()` function; must include a
  valid `learner` so the corresponding `predict_<learner>()` can be
  dispatched.

- newdata:

  A data frame with predictors for the subjects to cluster.

- times:

  Numeric vector of evaluation times (the common time grid on which
  survival curves are predicted and clustered).

- K:

  Optional integer number of clusters. If `NULL`, selected by silhouette
  width (see
  [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html)).

- K_max:

  Maximum `K` considered when `K` is `NULL`.

- distance:

  Distance type passed to
  [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html):
  `"L2"` or `"L1"`.

- seed:

  Optional integer seed forwarded to
  [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html).

- ...:

  Additional arguments forwarded to
  [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html) (e.g.
  `weights`, `enforce_monotone`, `smooth_median_width`,
  `standardize_cols`).

## Value

A list with components:

- fit:

  The fitted `unsurv` object (clusters, medoids, silhouette, ...).

- survmat:

  The numeric survival-probability matrix that was clustered.

- times:

  The evaluation time grid used.

- learner:

  The learner label (from `model$learner`).

## Details

Predicted survival is obtained from the model's `predict_<learner>()` on
`times`, coerced to a numeric matrix, and passed to
[`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html).
Requires the unsurv package.

## References

El Badisy, I. (2026). unsurv: Clustering Individualized Survival Curves.
Bioinformatics Advances, vbag218.
[doi:10.1093/bioadv/vbag218](https://doi.org/10.1093/bioadv/vbag218)

## See also

[`plot_survcluster()`](https://elbadisyimad.com/survalis/reference/plot_survcluster.md)

## Examples

``` r
if (requireNamespace("unsurv", quietly = TRUE) &&
    requireNamespace("cluster", quietly = TRUE)) {
  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
  time_grid <- seq(10, 300, length.out = 20)
  cl <- compute_survcluster(mod, newdata = veteran, times = time_grid, K = 2, seed = 1)
  table(cl$fit$clusters)
}
#> 
#>  1  2 
#> 67 70 
```
