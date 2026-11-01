# Fit a densemlp Survival Model

Fits a compact feedforward neural network for right-censored survival
data with the densemlp package (native C++ backend, no torch dependency)
and returns an object compatible with the `mlsurv_model` interface. The
IPCW integrated-Brier loss (`loss = "brier"`) is used by default so that
a full survival-probability curve is available; predictions are
interpolated from the network's discrete-time grid onto the requested
evaluation times.

## Usage

``` r
fit_densemlp(
  formula,
  data,
  loss = c("brier", "cox"),
  n_bins = 10,
  hidden_units = c(32, 16),
  epochs = 100,
  lr = 0.001,
  ...
)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the formula variables.

- loss:

  Survival loss: `"brier"` (default, discrete-time hazard grid, required
  for survival-curve prediction) or `"cox"` (a single linear risk score
  only).

- n_bins:

  Number of discrete-time bins for `loss = "brier"` (default `10`).

- hidden_units:

  Integer vector of hidden-layer widths (default `c(32, 16)`).

- epochs:

  Maximum training epochs (default `100`).

- lr:

  Adam learning rate (default `1e-3`).

- ...:

  Additional arguments passed to
  [`densemlp::densemlp()`](https://rdrr.io/pkg/densemlp/man/densemlp.html)
  (for example `batch_size`, `dropout`, `batch_norm`, `early_stopping`,
  `patience`, `ensemble`, `seed`, `ncores`).

## Value

An object of class `"mlsurv_model"` with elements `model`, `learner`
(`"densemlp"`), `formula`, `data`, `time`, `status`, and `loss`.

## See also

[`predict_densemlp`](https://elbadisyimad.com/survalis/reference/predict_densemlp.md),
[`densemlp`](https://rdrr.io/pkg/densemlp/man/densemlp.html),
[`fit`](https://elbadisyimad.com/survalis/reference/fit.md) for the verb
interface

## Examples

``` r
mod <- fit_densemlp(Surv(time, status) ~ age + karno + celltype, veteran,
                    n_bins = 8, epochs = 20, hidden_units = c(16, 8))
summary(mod)
#> 
#> ── densemlp summary ────────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: densemlp
#> Learner: densemlp
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
