#' Tidy data frames from verb result objects
#'
#' `as.data.frame()` methods that return the most useful rectangular view of
#' each result object, for downstream plotting or reporting:
#'
#' * `survalis_fit` -- a one-row model summary (`model`, `engine`, `formula`,
#'   `n`, `events`).
#' * `survmat` (from [predict()]) -- by default the wide prediction matrix with
#'   one `"t=<time>"` column per evaluation time; `format = "long"` returns one
#'   row per observation and time (`id`, `time`, `value`).
#' * `survalis_eval` -- the per-split metric values (`id`, `metric`, `value`).
#' * `survalis_tune` -- the scored hyperparameter grid.
#' * `survalis_compare` -- the model-by-metric leaderboard (`model`, `metric`,
#'   `mean`, `se`).
#' * `survalis_estimate` -- a one-row summary (`estimand`, `tau`, `estimate`,
#'   `conf_low`, `conf_high`, `model`).
#'
#' @param x A verb result object.
#' @param format For `survmat` only: `"wide"` (default) or `"long"`.
#' @param ... Ignored.
#' @return A data frame.
#' @name survalis-tidiers
NULL

#' @rdname survalis-tidiers
#' @export
as.data.frame.survalis_fit <- function(x, ...) {
  data.frame(model = x$model_id,
             engine = if (is.null(x$engine)) NA_character_ else as.character(x$engine),
             formula = paste(deparse(x$formula), collapse = " "),
             n = nrow(x$data),
             events = sum(.fit_status_vector(x)),
             stringsAsFactors = FALSE)
}

#' @rdname survalis-tidiers
#' @export
as.data.frame.survalis_eval <- function(x, ...) x$per_split

#' @rdname survalis-tidiers
#' @export
as.data.frame.survalis_tune <- function(x, ...) x$results

#' @rdname survalis-tidiers
#' @export
as.data.frame.survalis_compare <- function(x, ...) x$leaderboard

#' @rdname survalis-tidiers
#' @export
as.data.frame.survalis_estimate <- function(x, ...) summary(x)
