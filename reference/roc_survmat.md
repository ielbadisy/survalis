# Time-Dependent ROC Curve from a Survival-Probability Matrix

Computes the full cumulative/dynamic time-dependent ROC curve
(sensitivity vs. specificity across risk-score thresholds) at a fixed
horizon `t_star`, using the same inverse-probability-of-censoring
weighted (IPCW) case/control definition as
[`auc_survmat`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)
(Uno et al. 2007; Heagerty and Zheng 2005). Its trapezoidal AUC closely
tracks
[`auc_survmat`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)'s
pairwise-concordance AUC (small numerical differences arise from tie
handling between the step-function ROC curve and the exact pairwise
comparison).

## Usage

``` r
roc_survmat(object, predicted, t_star)

# S3 method for class 'roc_survmat'
plot(x, title, ...)

plot_roc(x, ...)
```

## Arguments

- object:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object of
  length \\n\\.

- predicted:

  An `n x k` matrix or data frame of survival probabilities with columns
  named `"t=<time>"`.

- t_star:

  Numeric time at which to evaluate the ROC curve.

- x:

  A `roc_survmat` object.

- title:

  Plot title. If missing, an automatically generated title is used
  (includes the AUC). Pass `NULL` to omit the title entirely.

- ...:

  Additional arguments passed to `plot.roc_survmat`.

## Value

An object of class `"roc_survmat"`: a list with `curve` (a data.frame
with `threshold`, `sensitivity`, `specificity`), `auc`, and `t_star`.

## Details

As in
[`auc_survmat`](https://elbadisyimad.com/survalis/reference/auc_survmat.md),
cases are subjects with an observed event strictly before `t_star`
(IPCW-weighted by \\1/\hat G(\tilde T_i)\\), and controls are subjects
known to survive beyond `t_star`. The risk score is \\1 - \hat S(t^\star
\mid x)\\. Sensitivity and specificity are swept over every observed
risk-score value \\c\\: \$\$ \widehat{\mathrm{Se}}(c) = \frac{\sum\_{i
\in \text{cases}} w_i\\ I(\text{risk}\_i \ge c)} {\sum\_{i \in
\text{cases}} w_i}, \qquad \widehat{\mathrm{Sp}}(c) = \frac{\sum\_{i \in
\text{controls}} I(\text{risk}\_i \< c)} {\|\text{controls}\|}. \$\$

## References

Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for
t-year survivors with censored regression models.
[doi:10.1198/016214507000000149](https://doi.org/10.1198/016214507000000149)

Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
curves.
[doi:10.1111/j.0006-341X.2005.030814.x](https://doi.org/10.1111/j.0006-341X.2005.030814.x)

## See also

[`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md),
`plot.roc_survmat()`

## Examples

``` r
y <- survival::Surv(time = veteran$time, event = veteran$status)
sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1,
            dimnames = list(NULL, "t=100"))
r <- roc_survmat(y, predicted = sp, t_star = 100)
r$auc
#> [1] 0.7962477
head(r$curve)
#>   threshold sensitivity specificity
#> 1 0.1173737   1.0000000   0.0000000
#> 2 0.1724406   1.0000000   0.0000000
#> 3 0.2109977   1.0000000   0.1320755
#> 4 0.2555142   0.9865336   0.1320755
#> 5 0.3057828   0.9090763   0.4339623
#> 6 0.3611437   0.8962462   0.4528302
```
