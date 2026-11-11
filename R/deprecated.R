#' Deprecated functions
#'
#' These names are kept for one release and will be removed in the next.
#' Use the glued names instead: [listlearners()], [listmetrics()],
#' [listmethods()] and [groupcv()].
#'
#' @param ... Passed to the new function.
#' @return The value of the new function.
#' @name survalis-deprecated
#' @keywords internal
NULL

#' @rdname survalis-deprecated
#' @export
list_survlearners <- function(...) {
  .Deprecated("listlearners", package = "survalis")
  listlearners(...)
}

#' @rdname survalis-deprecated
#' @export
list_metrics <- function(...) {
  .Deprecated("listmetrics", package = "survalis")
  listmetrics(...)
}

#' @rdname survalis-deprecated
#' @export
list_interpretability_methods <- function(...) {
  .Deprecated("listmethods", package = "survalis")
  listmethods(...)
}

#' @rdname survalis-deprecated
#' @export
group_cv <- function(...) {
  .Deprecated("groupcv", package = "survalis")
  groupcv(...)
}
