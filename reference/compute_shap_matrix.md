# Time-Integrated SHAP Values for a Set of Observations

Computes
[`compute_shap`](https://elbadisyimad.com/survalis/reference/compute_shap.md)
for every row of `newdata` and, instead of averaging across observations
(as the internal `compute_shap_mean()` helper does), keeps one
time-integrated SHAP value per (observation, feature) pair together with
each observation's raw feature value – the input needed for
[`plot_shap_beeswarm`](https://elbadisyimad.com/survalis/reference/plot_shap_beeswarm.md).

## Usage

``` r
compute_shap_matrix(
  model,
  newdata,
  baseline_data,
  times,
  sample.size = 100,
  features = NULL
)
```

## Arguments

- model:

  An `mlsurv_model` object from a `fit_*()` learner.

- newdata:

  A data.frame of observations to explain (one row per observation).

- baseline_data:

  A data.frame to sample background values from (typically the training
  data used to fit `model`).

- times:

  Numeric vector of evaluation times passed to
  [`compute_shap`](https://elbadisyimad.com/survalis/reference/compute_shap.md).

- sample.size:

  Integer, number of random feature orderings sampled per observation
  and time point (passed to
  [`compute_shap`](https://elbadisyimad.com/survalis/reference/compute_shap.md)).
  Default `100`.

- features:

  Optional character vector of covariate names to explain, passed
  through to
  [`compute_shap`](https://elbadisyimad.com/survalis/reference/compute_shap.md).
  If `NULL` (default), inferred as
  `setdiff(colnames(newdata), c(model$time, model$status))`.

## Value

A long data.frame with columns `observation`, `feature`, `phi`
(time-integrated SHAP value), and `raw_value` (the observation's value
for that feature).

## Details

For each observation \\x\\, per-time Shapley values \\\varphi_j(x, t)\\
(from
[`compute_shap`](https://elbadisyimad.com/survalis/reference/compute_shap.md))
are time-integrated via the trapezoidal rule over the requested time
grid, \$\$ \Phi_j(x) = \int_0^{T} \varphi_j(x, t)\\ dt \approx
\sum\_{k=2}^{m} \frac{\varphi_j(x, t_k) + \varphi_j(x, t\_{k-1})}{2}
(t_k - t\_{k-1}), \$\$ using the same integration already implemented by
`compute_shap(..., aggregate = TRUE, method = "integral")`. This
produces one summary contribution per (observation, feature) that
reflects impact over the *entire* prediction horizon in a single value,
which is what
[`plot_shap_beeswarm`](https://elbadisyimad.com/survalis/reference/plot_shap_beeswarm.md)
visualizes – as opposed to a per-timepoint SurvSHAP(t) panel.

## See also

[`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md),
[`plot_shap_beeswarm()`](https://elbadisyimad.com/survalis/reference/plot_shap_beeswarm.md)

## Examples

``` r
# \donttest{
mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
sm <- compute_shap_matrix(mod, newdata = veteran[1:15, ],
                          baseline_data = veteran, times = c(50, 100, 200),
                          sample.size = 5)
head(sm)
#>   observation  feature         phi raw_value
#> 1           1      trt  0.00000000         1
#> 2           1 celltype  0.00000000         1
#> 3           1    karno  1.42998911        60
#> 4           1 diagtime -0.00466043         7
#> 5           1      age  0.02071794        69
#> 6           1    prior  0.00000000         0
# }
```
