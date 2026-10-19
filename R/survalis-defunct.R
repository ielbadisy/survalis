#' Defunct functions in survalis
#'
#' These functions were removed from survalis and now only stop with a
#' pointer to their replacement.
#'
#' \tabular{ll}{
#'   \strong{Defunct} \tab \strong{Use instead} \cr
#'   `fit_survdnn()` \tab `fit(formula, data, model = "densemlp")` \cr
#'   `predict_survdnn()` \tab `predict()` on a `"densemlp"` fit \cr
#'   `tune_survdnn()` \tab `fit(formula, data, model = "densemlp", spec = ...)` \cr
#' }
#'
#' The `survdnn` learner was removed in survalis 1.7.0. It needed the
#' \pkg{torch} runtime, a large separate download that often failed to
#' install. The `densemlp` learner covers the same deep survival models with a
#' native C++ backend.
#'
#' @param ... Ignored.
#'
#' @name survalis-defunct
#' @keywords internal
NULL

.survdnn_defunct <- function(fun, new) {
  .Defunct(
    new = new, package = "survalis",
    msg = paste0(
      "`", fun, "()` was removed in survalis 1.7.0 together with the torch ",
      "dependency. Use the densemlp learner instead: ",
      "fit(formula, data, model = \"densemlp\")."
    )
  )
}

#' @rdname survalis-defunct
#' @export
fit_survdnn <- function(...) .survdnn_defunct("fit_survdnn", "fit_densemlp")

#' @rdname survalis-defunct
#' @export
predict_survdnn <- function(...) .survdnn_defunct("predict_survdnn", "predict_densemlp")

#' @rdname survalis-defunct
#' @export
tune_survdnn <- function(...) .survdnn_defunct("tune_survdnn", "fit_densemlp")
