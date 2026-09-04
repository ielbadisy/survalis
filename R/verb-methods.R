#' Tidy data frames from verb result objects
#'
#' `as.data.frame()` methods that return the most useful rectangular view of
#' each result object, for downstream plotting or reporting:
#'
#' * `survalis_eval` -- the per-split metric values (`id`, `metric`, `value`).
#' * `survalis_tune` -- the scored hyperparameter grid.
#' * `survalis_compare` -- the model-by-metric leaderboard (`model`, `metric`,
#'   `mean`, `se`).
#' * `survalis_estimate` -- a one-row summary (`estimand`, `tau`, `estimate`,
#'   `conf.low`, `conf.high`, `model`).
#'
#' @param x A verb result object.
#' @param ... Ignored.
#' @return A data frame.
#' @name survalis-tidiers
NULL

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
