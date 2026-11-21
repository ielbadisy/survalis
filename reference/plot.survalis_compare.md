# Plot a model comparison

One panel per metric, with the per-split values of each model as a box
plot and points. Each panel has its own value axis, and its strip says
whether higher or lower is better. Models are ordered by the first
metric, best on top.

## Usage

``` r
# S3 method for class 'survalis_compare'
plot(x, title, scales = "free_x", ...)
```

## Arguments

- x:

  A `survalis_compare` object.

- title:

  Plot title; `NULL` for none.

- scales:

  Passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html):
  `"free_x"` (default) gives each metric its own value axis, `"fixed"`
  one shared axis.

- ...:

  Unused.

## Value

A ggplot.
