# Fit a Survival Support Vector Machine

Fits a kernel survival support vector regression on the log-time scale
and returns an `mlsurv_model`. The engine is implemented natively in
survalis.

## Usage

``` r
fit_survsvm(
  formula,
  data,
  type = "regression",
  gamma.mu = 0.1,
  opt.meth = NULL,
  kernel = c("add_kernel", "lin_kernel", "rbf_kernel", "poly_kernel"),
  diff.meth = NULL,
  kernel.pars = NULL,
  resid.folds = 5
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- type:

  Only `"regression"` is supported: the ranking types of survivalsvm
  give a risk score, not a survival time, and therefore no survival
  curve.

- gamma.mu:

  Error weight \\\gamma \> 0\\; larger values follow the data more
  closely. Default `0.1`.

- opt.meth, diff.meth:

  Ignored. Kept so that older calls keep working; the native solver
  needs no optimizer choice.

- kernel:

  Kernel: `"add_kernel"` (default), `"lin_kernel"`, `"rbf_kernel"` or
  `"poly_kernel"`.

- kernel.pars:

  Optional kernel parameters, see the Kernels section.

- resid.folds:

  Number of folds for the cross-fitted residuals (default `5`); `1` uses
  in-sample residuals. Folds interleave the ordered survival times, so
  fitting does not touch the random number generator.

## Value

An object of class `mlsurv_model`, a named list with elements `model`
(the fitted engine), `learner` (`"survsvm"`), `formula`, `data`, `time`
and `status`; the engine name `"survalis"` is stored in
`attr(, "engine")`.

## Method

The model is an accelerated failure time model \\\log T = f(x) +
\varepsilon\\ with \\f(x) = \sum_i \beta_i k(x, x_i) + b\\ in the
reproducing kernel Hilbert space of the kernel \\k\\. With \\y_i = \log
T_i\\, \\f\\ minimizes \$\$\frac12 \\f\\\_k^2 + \frac{\gamma}{2}
\sum\_{i:\\\delta_i = 1} (y_i - f(x_i))^2 + \frac{\gamma}{2}
\sum\_{i:\\\delta_i = 0} \max(0,\\ y_i - f(x_i))^2 ,\$\$ the survival
support vector regression of Shivaswamy et al. (2007) and Van Belle et
al. (2011) with squared losses: an event time is fitted from both sides,
a censored time is only penalized when the prediction falls below it.
The problem is strictly convex and piecewise quadratic, and is solved
exactly by an active-set Newton method (Chapelle, 2007). Given the
active set \\A\\ (events, and censored points predicted too early), the
optimum solves the least-squares SVM system \$\$\begin{pmatrix} 0 &
1^\top \\ 1 & K\_{AA} + I/\gamma \end{pmatrix} \begin{pmatrix} b \\
\beta_A \end{pmatrix} = \begin{pmatrix} 0 \\ y_A \end{pmatrix}, \qquad
\beta\_{A^c} = 0,\$\$ by one Cholesky factorization; \\A\\ is updated,
with a line search, until it no longer changes. No quadratic programming
solver is involved. Censored points outside \\A\\ have \\\beta_i = 0\\
and are not stored.

Survival curves follow from the AFT form: \\S(t \mid x) =
S\_\varepsilon(\log t - f(x))\\, where \\S\_\varepsilon\\ is the
Kaplan-Meier estimate of the cross-fitted residuals \\y_i -
f\_{-k(i)}(x_i)\\: each residual comes from a model fitted without the
observation's fold, so the curves are not made overconfident by the
in-sample fit. See
[`predict_survsvm`](https://elbadisyimad.com/survalis/reference/predict_survsvm.md)
for alternatives.

## Kernels

- `"add_kernel"`:

  Clinical (additive) kernel of Daemen et al. (2012): the average over
  predictors of \\\max(0, 1 - \|x_j - z_j\| / r_j)\\ for numeric and
  ordered predictors (\\r_j\\ the training range) and \\1\\x_j = z_j\\\\
  for unordered factors. Scale free, positive semi-definite, and a
  factor counts as one predictor. Constant predictors add a constant.

- `"lin_kernel"`:

  \\x^\top z / p\\ on standardized columns (factors in treatment
  coding).

- `"rbf_kernel"`:

  \\\exp(-\\x - z\\^2 / (2 \sigma^2))\\ on standardized columns;
  `kernel.pars` is \\\sigma^2\\, default the number of columns \\p\\.

- `"poly_kernel"`:

  \\(x^\top z / p + c)^d\\ on standardized columns; `kernel.pars` is
  \\(c, d)\\, default \\(1, 2)\\.

## References

Shivaswamy PK, Chu W, Jansche M (2007). A support vector approach to
censored targets. *Seventh IEEE International Conference on Data
Mining*, 655-660.
[doi:10.1109/ICDM.2007.93](https://doi.org/10.1109/ICDM.2007.93)

Van Belle V, Pelckmans K, Van Huffel S, Suykens JAK (2011). Support
vector methods for survival analysis: a comparison between ranking and
regression approaches. *Artificial Intelligence in Medicine*, 53(2),
107-118.
[doi:10.1016/j.artmed.2011.06.006](https://doi.org/10.1016/j.artmed.2011.06.006)

Daemen A, Timmerman D, Van den Bosch T, et al. (2012). Improved modeling
of clinical data with kernel methods. *Artificial Intelligence in
Medicine*, 54(2), 103-114.
[doi:10.1016/j.artmed.2011.11.001](https://doi.org/10.1016/j.artmed.2011.11.001)

Chapelle O (2007). Training a support vector machine in the primal.
*Neural Computation*, 19(5), 1155-1178.
[doi:10.1162/neco.2007.19.5.1155](https://doi.org/10.1162/neco.2007.19.5.1155)

## See also

[`predict_survsvm()`](https://elbadisyimad.com/survalis/reference/predict_survsvm.md),
[`tune_survsvm()`](https://elbadisyimad.com/survalis/reference/tune_survsvm.md)

## Examples

``` r
mod_svm <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
                       data = veteran, kernel = "add_kernel")
predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = c(100, 300, 500))
#>       t=100      t=300      t=500
#> 1 0.4911426 0.11209330 0.04670554
#> 2 0.5367487 0.14887392 0.06538776
#> 3 0.4450980 0.09341109 0.04670554
#> 4 0.4988167 0.11209330 0.04670554
#> 5 0.5367487 0.14887392 0.06538776
```
