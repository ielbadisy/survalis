# FDR-Based Variable Screening via Univariate Cox Models

Screens candidate covariates for a right-censored outcome by fitting a
univariate Cox model \\\lambda(t \mid X_j) = \lambda_0(t)\exp(\beta_j
X_j)\\ per covariate, Wald-testing \\H_0: \beta_j = 0\\, and controlling
the false discovery rate (FDR) across all tests via the
Benjamini-Hochberg (1995) step-up procedure.

## Usage

``` r
screen_fdr(formula, data, alpha = 0.1, minscreen = 5)
```

## Arguments

- formula:

  A model formula with a `Surv(time, status)` outcome on the left-hand
  side and candidate screening covariates on the right-hand side (as
  used throughout survalis, e.g.
  [`fit_coxph`](https://elbadisyimad.com/survalis/reference/fit_coxph.md)).

- data:

  A data.frame containing the variables in `formula`.

- alpha:

  Target FDR level in (0, 1). Covariates with BH-adjusted p-value
  (q-value) `<= alpha` are retained. Default `0.1`.

- minscreen:

  Minimum number of covariates to retain. If fewer than `minscreen`
  covariates pass the FDR threshold, the smallest-p-value covariates are
  added back until `minscreen` is reached (subject to the number of
  candidate covariates available). Default `5`.

## Value

A character vector of selected covariate names, with an attribute
`"screen_table"` holding the full `data.frame` of `feature`, `p_value`,
`q_value`, and `selected`.

## Details

For each candidate covariate \\X_j\\, survalis fits
`coxph(Surv(time, status) ~ X_j, data = data)` and extracts the Wald
p-value for \\\beta_j\\. Raw p-values \\p_1, \dots, p_p\\ are ordered
\\p\_{(1)} \le \cdots \le p\_{(p)}\\ and the BH-adjusted q-value for the
\\k\\-th smallest p-value is \$\$q\_{(k)} = \min\_{k \le l \le p}
\left\\ \frac{p \cdot p\_{(l)}}{l} \right\\,\$\$ guaranteeing
\\\mathrm{E}\[\text{FDR}\] \le \alpha\\ under the retained set \\\\j :
q_j \le \alpha\\\\ (Benjamini and Hochberg 1995).

## References

Benjamini Y, Hochberg Y (1995). "Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing." *Journal of the
Royal Statistical Society: Series B*, 57(1), 289-300.

## Examples

``` r
selected <- screen_fdr(Surv(time, status) ~ age + karno + diagtime + prior,
                       data = veteran, alpha = 0.2, minscreen = 2)
selected
#> [1] "diagtime" "karno"   
#> attr(,"screen_table")
#>    feature      p_value      q_value selected
#> 1      age 4.329948e-01 4.759247e-01    FALSE
#> 2 diagtime 3.107793e-01 4.759247e-01     TRUE
#> 3    karno 4.513209e-11 1.805284e-10     TRUE
#> 4    prior 4.759247e-01 4.759247e-01    FALSE
attr(selected, "screen_table")
#>    feature      p_value      q_value selected
#> 1      age 4.329948e-01 4.759247e-01    FALSE
#> 2 diagtime 3.107793e-01 4.759247e-01     TRUE
#> 3    karno 4.513209e-11 1.805284e-10     TRUE
#> 4    prior 4.759247e-01 4.759247e-01    FALSE
```
