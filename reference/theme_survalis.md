# A Consistent ggplot2 Theme for survalis Figures

A clean, minimal theme shared across survalis plotting functions
([`plot_pdp`](https://elbadisyimad.com/survalis/reference/plot_pdp.md),
[`plot_benchmark`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md),
[`cv_plot`](https://elbadisyimad.com/survalis/reference/cv_plot.md), and
others), so figures produced by the package look like one consistent
system rather than each using an independent ad-hoc style.

## Usage

``` r
theme_survalis(base_size = 13, legend.position = "bottom")
```

## Arguments

- base_size:

  Base font size (default `13`).

- legend.position:

  Legend position passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  (default `"bottom"`).

## Value

A ggplot2 theme object, to be added to a ggplot with `+`.

## See also

[`scale_color_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md),
[`scale_fill_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)

## Examples

``` r
p <- ggplot2::ggplot(veteran, ggplot2::aes(age, karno)) + ggplot2::geom_point()
p + theme_survalis()

```
