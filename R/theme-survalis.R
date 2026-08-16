.survalis_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666"
)

.survalis_pal <- function(n) {
  if (n <= length(.survalis_palette)) {
    .survalis_palette[seq_len(n)]
  } else {
    grDevices::colorRampPalette(.survalis_palette)(n)
  }
}

#' A Consistent ggplot2 Theme for survalis Figures
#'
#' A clean, minimal theme shared across \pkg{survalis} plotting functions
#' (\code{\link{plot_pdp}}, \code{\link{plot_benchmark}}, \code{\link{cv_plot}},
#' and others), so figures produced by the package look like one consistent
#' system rather than each using an independent ad-hoc style.
#'
#' @param base_size Base font size (default \code{13}).
#' @param legend.position Legend position passed to \code{ggplot2::theme()}
#'   (default \code{"bottom"}).
#'
#' @return A \pkg{ggplot2} theme object, to be added to a ggplot with \code{+}.
#'
#' @examples
#' p <- ggplot2::ggplot(veteran, ggplot2::aes(age, karno)) + ggplot2::geom_point()
#' p + theme_survalis()
#'
#' @seealso [scale_color_survalis()], [scale_fill_survalis()]
#' @export
theme_survalis <- function(base_size = 13, legend.position = "bottom") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      axis.title = ggplot2::element_text(face = "plain"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = legend.position,
      legend.title = ggplot2::element_text(face = "plain")
    )
}

#' Discrete Color/Fill Scales Matching \code{\link{theme_survalis}}
#'
#' Colorblind-friendly discrete color and fill scales (based on ColorBrewer
#' "Dark2", recycled/interpolated beyond 8 levels) for consistent group
#' colors across \pkg{survalis} figures.
#'
#' @param ... Additional arguments passed to \code{ggplot2::discrete_scale()}.
#'
#' @return A \pkg{ggplot2} discrete scale object.
#'
#' @examples
#' p <- ggplot2::ggplot(veteran, ggplot2::aes(age, karno, color = celltype)) +
#'   ggplot2::geom_point()
#' p + scale_color_survalis()
#'
#' @name survalis-scales
#' @export
scale_color_survalis <- function(...) {
  ggplot2::discrete_scale("colour", palette = .survalis_pal, ...)
}

#' @rdname survalis-scales
#' @export
scale_fill_survalis <- function(...) {
  ggplot2::discrete_scale("fill", palette = .survalis_pal, ...)
}
