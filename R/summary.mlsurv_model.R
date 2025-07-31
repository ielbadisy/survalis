#' Summarize a Survival Learner Model (`mlsurv_model`)
#'
#' Provides a structured overview of a fitted survival model, including formula,
#' engine, and data characteristics. Styled like `summary.survdnn()` and supports
#' all learners wrapped using the `mlsurv_model` format.
#'
#' @param object An object of class `"mlsurv_model"` returned by a `fit_*` function.
#' @param ... Currently ignored (for future compatibility).
#'
#' @return Invisibly returns an object of class `"summary.mlsurv_model"`.
#' @export
#'
#' @examples
#' data(veteran)
#' mod <- fit_cox(Surv(time, status) ~ age + trt + celltype, data = veteran)
#' summary(mod)
summary.mlsurv_model <- function(object, ...) {
  stopifnot(inherits(object, "mlsurv_model"))

  mf <- model.frame(object$formula, object$data)
  xmat <- model.matrix(delete.response(terms(object$formula)), mf)
  xnames <- setdiff(colnames(xmat), "(Intercept)")
  y <- model.response(mf)
  time <- y[, "time"]
  status <- y[, "status"]

  out <- list(
    learner = object$learner,
    engine = attr(object, "engine"),
    formula = object$formula,
    data_summary = list(
      observations = nrow(object$data),
      predictors = xnames,
      time_range = range(time, na.rm = TRUE),
      event_rate = mean(status, na.rm = TRUE)
    )
  )
  class(out) <- "summary.mlsurv_model"

  # styled output
  cli::cli_h1(paste("Summary of", out$learner, "survival model"))

  cli::cli_text("Formula:")
  cli::cli_code(deparse(out$formula))

  cli::cli_text("Engine: {.strong {out$engine}}")
  cli::cli_text("Learner: {.strong {out$learner}}")

  cli::cli_text("\nData summary:")
  cli::cli_text("- Observations: {.val {out$data_summary$observations}}")
  cli::cli_text("- Predictors: {.val {paste(out$data_summary$predictors, collapse = ', ')}}")
  cli::cli_text("- Time range: [{out$data_summary$time_range[1]}, {out$data_summary$time_range[2]}]")
  cli::cli_text("- Event rate: {.val {sprintf('%.1f%%', 100 * out$data_summary$event_rate)}}")

  invisible(out)
}
