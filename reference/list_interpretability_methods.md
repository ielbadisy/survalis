# List interpretability methods available in survalis

Returns a table mapping each `compute_*` function to its paired `plot_*`
helper (if any). Methods without a plot helper show `NA`.

## Usage

``` r
list_interpretability_methods()
```

## Value

A data.frame with columns `compute`, `plot`, `has_compute`, and
`has_plot`.

## Examples

``` r
list_interpretability_methods()
#>                                       compute                plot has_compute
#> compute_shap                     compute_shap           plot_shap        TRUE
#> compute_pdp                       compute_pdp            plot_pdp        TRUE
#> compute_ale                       compute_ale            plot_ale        TRUE
#> compute_surrogate           compute_surrogate      plot_surrogate        TRUE
#> compute_tree_surrogate compute_tree_surrogate plot_tree_surrogate        TRUE
#> compute_varimp                 compute_varimp         plot_varimp        TRUE
#> compute_interactions     compute_interactions   plot_interactions        TRUE
#> compute_counterfactual compute_counterfactual plot_counterfactual        TRUE
#> compute_shap_matrix       compute_shap_matrix  plot_shap_beeswarm        TRUE
#> compute_survcluster       compute_survcluster    plot_survcluster        TRUE
#>                        has_plot
#> compute_shap               TRUE
#> compute_pdp                TRUE
#> compute_ale                TRUE
#> compute_surrogate          TRUE
#> compute_tree_surrogate     TRUE
#> compute_varimp             TRUE
#> compute_interactions       TRUE
#> compute_counterfactual     TRUE
#> compute_shap_matrix        TRUE
#> compute_survcluster        TRUE
subset(list_interpretability_methods(), is.na(plot))  # methods without plot
#> [1] compute     plot        has_compute has_plot   
#> <0 rows> (or 0-length row.names)
```
