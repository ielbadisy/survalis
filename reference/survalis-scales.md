# Discrete Color/Fill Scales Matching [`theme_survalis`](https://elbadisyimad.com/survalis/reference/theme_survalis.md)

Colorblind-friendly discrete color and fill scales (based on ColorBrewer
"Dark2", recycled/interpolated beyond 8 levels) for consistent group
colors across survalis figures.

## Usage

``` r
scale_color_survalis(...)

scale_fill_survalis(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

## Value

A ggplot2 discrete scale object.

## Examples

``` r
p <- ggplot2::ggplot(veteran, ggplot2::aes(age, karno, color = celltype)) +
  ggplot2::geom_point()
p + scale_color_survalis()

```
