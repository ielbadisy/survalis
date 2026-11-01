# Decision Curve Analysis for Right-Censored Survival Outcomes

Computes net benefit across a range of risk thresholds at a fixed
horizon `t_star` for a survival model's predicted risk, alongside the
"treat all" and "treat none" reference strategies (Vickers and Elkin
2006), adapted to right-censored outcomes via the same
inverse-probability-of- censoring weighting (IPCW) estimator of \\\hat
G(t)\\ used by
[`brier`](https://elbadisyimad.com/survalis/reference/brier.md) and
[`auc_survmat`](https://elbadisyimad.com/survalis/reference/auc_survmat.md).

## Usage

``` r
dca_survmat(object, predicted, t_star, thresholds = seq(0.01, 0.99, by = 0.01))

# S3 method for class 'dca_survmat'
plot(x, title, ...)

plot_dca(x, ...)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object of
  length \\n\\.

- predicted:

  An `n x k` matrix or data frame of survival probabilities with columns
  named `"t=<time>"`.

- t_star:

  Numeric horizon at which risk is evaluated.

- thresholds:

  Numeric vector of risk thresholds in (0, 1). Default
  `seq(0.01, 0.99, by = 0.01)`.

- x:

  A `dca_survmat` object.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely.

- ...:

  Additional arguments passed to `plot.dca_survmat`.

## Value

An object of class `"dca_survmat"`: a list with `curve` (a data.frame
with `threshold`, `model`, `treat_all`, `treat_none`) and `t_star`.

## Details

The predicted risk at the horizon is \\\pi_i = 1 - \hat S(t^\star \mid
X_i)\\. At threshold `pt`, subjects with \\\pi_i \ge pt\\ are classified
high-risk. Using the IPCW identity for the Kaplan-Meier estimator (Graf
et al. 1999), the true- and false-positive fractions are estimated as
\$\$ \widehat{\mathrm{TPF}}(pt) = \frac{1}{n}\sum\_{i=1}^n I(\pi_i \ge
pt)\\ \frac{\Delta_i\\ I(\tilde T_i \le t^\star)}{\hat G(\tilde T_i)},
\qquad \widehat{\mathrm{FPF}}(pt) = \frac{1}{n}\sum\_{i=1}^n I(\pi_i \ge
pt)\\ \frac{I(\tilde T_i \> t^\star)}{\hat G(t^\star)}, \$\$ and the net
benefit is \$\$ \widehat{\mathrm{NB}}(pt) = \widehat{\mathrm{TPF}}(pt) -
\widehat{\mathrm{FPF}}(pt) \cdot \frac{pt}{1-pt}. \$\$ "Treat all"
(classify everyone high-risk) reduces to the marginal IPCW-KM identity
\\1 - \hat S\_{\mathrm{KM}}(t^\star)\\ for the TPF term and \\\hat
S\_{\mathrm{KM}}(t^\star)\\ for the FPF term. "Treat none" has net
benefit \\0\\ by definition.

## References

Vickers AJ, Elkin EB (2006). "Decision Curve Analysis: A Novel Method
for Evaluating Prediction Models." *Medical Decision Making*, 26(6),
565-574.

Vickers AJ, Van Calster B, Steyerberg EW (2019). "A Simple, Step-by-Step
Guide to Interpreting Decision Curve Analysis." *Diagnostic and
Prognostic Research*, 3, 18.

Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). "Assessment and
Comparison of Prognostic Classification Schemes for Survival Data."
*Statistics in Medicine*, 18(17-18), 2529-2545.

## See also

`plot.dca_survmat()`,
[`brier()`](https://elbadisyimad.com/survalis/reference/brier.md),
[`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1,
            dimnames = list(NULL, "t=100"))
d <- dca_survmat(y, predicted = sp, t_star = 100)
head(d$curve)
#>   threshold     model  treat_all treat_none
#> 1      0.01 0.5781586 0.06430204          0
#> 2      0.02 0.5738541 0.05475410          0
#> 3      0.03 0.5694608 0.04500930          0
#> 4      0.04 0.5649760 0.03506148          0
#> 5      0.05 0.5603968 0.02490423          0
#> 6      0.06 0.5557202 0.01453087          0
```
