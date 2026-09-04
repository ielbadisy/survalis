#' The survmat prediction contract
#'
#' Internal constructor and validator for `survmat`, the single prediction
#' container used by the verb-level [predict()] method. A `survmat` is an
#' `n x length(times)` numeric matrix carrying its evaluation grid in the
#' `"times"` attribute (column names `"t=<time>"` are kept for display and for
#' backward compatibility with the `*_survmat` metric helpers) and the predicted
#' quantity in the `"quantity"` attribute.
#'
#' `new_survmat()` builds the object without checks; `validate_survmat()` enforces
#' the contract; `as_survmat()` is the checked constructor used internally after a
#' learner's raw `predict_*()` output has been finalised.
#'
#' @param x A numeric matrix (or data frame / vector coercible to one) of
#'   predictions, one row per observation and one column per entry of `times`.
#' @param times Numeric vector, strictly increasing and non-negative, of length
#'   `ncol(x)`.
#' @param quantity One of `"survival"`, `"risk"`, `"chf"`, `"hazard"`. Controls
#'   printing and the `[0, 1]` range check (applied to `"survival"` and
#'   `"risk"` only).
#'
#' @return A `survmat`: a numeric matrix with class `c("survmat", "matrix",
#'   "array")`, attribute `"times"`, and attribute `"quantity"`.
#' @keywords internal
#' @name survmat
NULL

#' @rdname survmat
#' @keywords internal
new_survmat <- function(x, times, quantity = "survival") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x) && !is.list(x)) x <- matrix(x, nrow = 1L)
  storage.mode(x) <- "double"
  times <- as.numeric(times)
  colnames(x) <- paste0("t=", times)
  rownames(x) <- NULL
  attr(x, "times") <- times
  attr(x, "quantity") <- quantity
  class(x) <- c("survmat", "matrix", "array")
  x
}

#' @rdname survmat
#' @keywords internal
validate_survmat <- function(x) {
  if (!inherits(x, "survmat") || !is.matrix(x) || !is.numeric(x)) {
    stop("`x` is not a valid survmat (expected a numeric matrix).", call. = FALSE)
  }
  times <- attr(x, "times")
  quantity <- attr(x, "quantity")
  if (is.null(times) || length(times) != ncol(x)) {
    stop("survmat `times` attribute must have length ncol(x).", call. = FALSE)
  }
  if (any(!is.finite(times)) || any(times < 0)) {
    stop("survmat `times` must be finite and non-negative.", call. = FALSE)
  }
  if (length(times) > 1L && is.unsorted(times, strictly = TRUE)) {
    stop("survmat `times` must be strictly increasing.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("survmat contains missing values; predictions must be complete ",
         "(check `newdata` for rows with NA in model variables).", call. = FALSE)
  }
  if (!is.null(quantity) && quantity %in% c("survival", "risk")) {
    rng <- range(x)
    if (rng[1] < -1e-8 || rng[2] > 1 + 1e-8) {
      stop("survmat with quantity '", quantity, "' has values outside [0, 1].",
           call. = FALSE)
    }
  }
  invisible(x)
}

#' @rdname survmat
#' @keywords internal
as_survmat <- function(x, times, quantity = "survival") {
  validate_survmat(new_survmat(x, times, quantity))
}

#' @keywords internal
is_survmat <- function(x) inherits(x, "survmat")

#' Evaluation grid of a survmat
#'
#' Returns the numeric time grid of a `survmat`, reading the `"times"` attribute
#' and falling back to parsing `"t=<time>"` column names for objects produced by
#' the legacy `predict_*()` path.
#'
#' @param x A `survmat` or a legacy survival-probability matrix / data frame.
#' @return A numeric vector.
#' @keywords internal
survmat_times <- function(x) {
  tt <- attr(x, "times")
  if (!is.null(tt)) return(as.numeric(tt))
  .infer_survmat_times(x)
}

#' @export
print.survmat <- function(x, n = 6L, ...) {
  q <- attr(x, "quantity")
  if (is.null(q)) q <- "survival"
  tt <- attr(x, "times")
  cat(sprintf("<survmat> %s  |  %d obs x %d times  [%s, %s]\n",
              q, nrow(x), ncol(x),
              format(min(tt)), format(max(tt))))
  body <- unclass(x)
  attr(body, "times") <- NULL
  attr(body, "quantity") <- NULL
  print(utils::head(body, n))
  if (nrow(x) > n) cat(sprintf("... %d more rows\n", nrow(x) - n))
  invisible(x)
}

#' @export
as.data.frame.survmat <- function(x, ...) {
  body <- unclass(x)
  attr(body, "times") <- NULL
  attr(body, "quantity") <- NULL
  as.data.frame(body, ...)
}

#' @export
as.matrix.survmat <- function(x, ...) {
  body <- unclass(x)
  attr(body, "times") <- NULL
  attr(body, "quantity") <- NULL
  body
}
