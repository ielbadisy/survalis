#' Summarize an `mlsurv_model`
#'
#' Produces a compact, human‑readable overview of a fitted survival learner that
#' follows the `mlsurv_model` contract (e.g., objects returned by `fit_*()` in
#' this toolkit). The summary prints the learner id, engine, original formula,
#' and basic data characteristics (sample size, predictor names, time range,
#' event rate). A structured list is returned invisibly for programmatic use.
#'
#' @param object An object of class `"mlsurv_model"` produced by a `fit_*()` function.
#' @param ... Ignored; included for S3 signature compatibility.
#'
#' @return Invisibly returns a list of class `"summary.mlsurv_model"` containing:
#' \describe{
#'   \item{learner}{Character id of the learner (e.g., `"ranger"`, `"coxph"`).}
#'   \item{engine}{Underlying package/engine used to fit the model.}
#'   \item{formula}{The original survival formula.}
#'   \item{data_summary}{List with \code{observations}, \code{predictors},
#'         \code{time_range}, and \code{event_rate}.}
#' }
#'
#' @details
#' This method relies on the presence of the fields \code{learner}, \code{formula},
#' and \code{data} stored in the fitted object (the standard `mlsurv_model`
#' contract). Output printing uses \pkg{cli} if available.
#'
#' @examples
#' mod <- fit_coxph(Surv(time, status) ~ age + trt + celltype, data = veteran)
#' summary(mod)
#'
#' s <- summary(mod)  # capture the structured result invisibly
#' str(s)
#'
#' @seealso [fit_coxph()], other `fit_*()` learners returning `mlsurv_model` objects.
#'
#' @export
#' @method summary mlsurv_model

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

  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_h1(paste(out$learner, "summary"))
    cli::cli_text("Formula:")
    cli::cli_code(deparse(out$formula))
    cli::cli_text("Engine: {.strong {out$engine}}")
    cli::cli_text("Learner: {.strong {out$learner}}")
    cli::cli_text("\nData summary:")
    cli::cli_text("- Observations: {.val {out$data_summary$observations}}")
    cli::cli_text("- Predictors: {.val {paste(out$data_summary$predictors, collapse = ', ')}}")
    cli::cli_text("- Time range: [{out$data_summary$time_range[1]}, {out$data_summary$time_range[2]}]")
    cli::cli_text("- Event rate: {.val {sprintf('%.1f%%', 100 * out$data_summary$event_rate)}}")
  } else {
    cat(paste0("\n-- ", out$learner, " summary --\n"))
    cat("Formula:\n")
    cat(paste0(deparse(out$formula), "\n"))
    cat("Engine:", out$engine, "\n")
    cat("Learner:", out$learner, "\n")
    cat("\nData summary:\n")
    cat("- Observations:", out$data_summary$observations, "\n")
    cat("- Predictors:", paste(out$data_summary$predictors, collapse = ", "), "\n")
    cat("- Time range: [", out$data_summary$time_range[1], ", ", out$data_summary$time_range[2], "]\n", sep = "")
    cat("- Event rate:", sprintf("%.1f%%", 100 * out$data_summary$event_rate), "\n")
  }

  invisible(out)
}
