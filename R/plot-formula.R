#' Plot a formula
#'
#' Dispatches on the left-hand side of the formula. A survival outcome
#' (`Surv(time, status) ~ group`) draws a styled Kaplan-Meier curve via
#' \code{\link{plot_survcurve}}; any other formula falls through to the
#' standard \code{graphics::plot.formula()} scatterplot method, so ordinary
#' \code{plot(y ~ x, data = df)} calls from other code are unaffected.
#'
#' @param formula A formula. \code{Surv(...) ~ ...} draws a KM curve via
#'   \code{\link{plot_survcurve}}; anything else is forwarded to
#'   \code{graphics::plot.formula()}.
#' @param data A data frame. Required for a survival formula.
#' @param ... For a survival formula, additional arguments to
#'   \code{\link{plot_survcurve}} (\code{conf.int}, \code{risk.table},
#'   \code{pval}, \code{title}, ...); otherwise forwarded to
#'   \code{graphics::plot.formula()}.
#'
#' @return For a survival formula, the \pkg{ggplot2} object returned by
#'   \code{\link{plot_survcurve}}. Otherwise, whatever
#'   \code{graphics::plot.formula()} returns (a base-graphics plot, drawn as a
#'   side effect).
#'
#' @seealso \code{\link{plot_survcurve}}
#'
#' @examples
#' plot(Surv(time, status) ~ trt, data = veteran)
#' plot(mpg ~ wt, data = mtcars)
#'
#' @method plot formula
#' @export
plot.formula <- function(formula, data = NULL, ...) {
  if (.is_surv_lhs(formula)) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("plot() on a Surv(...) formula requires a data frame `data`.", call. = FALSE)
    }
    return(plot_survcurve(formula, data = data, ...))
  }

  # Not a survival formula: defer to the base graphics method. Looked up by
  # name inside the graphics namespace (not via the S3 dispatch table) so this
  # always resolves to the original implementation, regardless of what other
  # packages have registered for plot.formula.
  base_plot_formula <- utils::getS3method("plot", "formula", envir = asNamespace("graphics"))
  if (is.null(data)) {
    base_plot_formula(formula, ...)
  } else {
    base_plot_formula(formula, data = data, ...)
  }
}

#' Does a formula have a Surv(...) left-hand side?
#'
#' Purely syntactic check (no `data` needed): mirrors the detection used by
#' \code{\link{.parse_surv_formula}} but does not require the formula to be
#' resolvable against a data frame.
#'
#' @param formula A formula.
#' @return A single logical.
#' @keywords internal
.is_surv_lhs <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) < 3L) return(FALSE)
  lhs <- formula[[2]]
  if (!is.call(lhs)) return(FALSE)
  head <- lhs[[1]]
  identical(head, as.name("Surv")) ||
    (is.call(head) && as.character(head[[1]]) %in% c("::", ":::") &&
       identical(head[[3]], as.name("Surv")))
}
