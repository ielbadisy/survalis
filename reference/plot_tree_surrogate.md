# Plot Tree-Based Surrogate Models or Feature Importances

Visualizes the results of a `tree_surrogate` object returned by
[`compute_tree_surrogate()`](https://elbadisyimad.com/survalis/reference/compute_tree_surrogate.md).
Can display either the fitted surrogate trees or aggregated feature
importance based on split counts.

## Usage

``` r
plot_tree_surrogate(
  tree_surrogate,
  type = c("tree", "importance"),
  top_n = 10,
  title
)
```

## Arguments

- tree_surrogate:

  An object of class `"tree_surrogate"`.

- type:

  Character string indicating the type of plot:

  - `"tree"`: displays the surrogate decision tree(s) for each time.

  - `"importance"`: displays a bar plot of top features ranked by split
    count.

- top_n:

  Integer, the number of top features to display in the importance plot.
  Ignored if `type = "tree"`.

- title:

  Plot title for `type = "importance"`. If missing, an automatically
  generated title is used. Pass `NULL` to omit the title entirely (e.g.,
  for journals requiring caption-only figures). Ignored for
  `type = "tree"`, which titles each per-time tree panel individually.

## Value

A plot object (for `"importance"`) or printed tree diagrams (for
`"tree"`).

## Details

- If `type = "tree"`, a separate tree diagram is produced for each
  evaluation time.

- If `type = "importance"`, feature split counts are summed across all
  times and plotted as a bar chart.

Requires the `partykit` package for tree plotting and `ggplot2` for
importance plotting.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_tree_surrogate(tree_ranger, type = "tree")
plot_tree_surrogate(tree_ranger, type = "importance", top_n = 5)
} # }
```
