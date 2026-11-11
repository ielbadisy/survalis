#' Fit a survival learner
#'
#' The entry point of the survalis verb interface. `fit()` selects a learner by
#' id, passes engine arguments through `spec`, and returns a `survalis_fit`
#' object that the other verbs ([predict()], [evaluate()], [tune()], [interpret()],
#' [estimate()]) consume.
#'
#' @param formula A survival formula with a `Surv()` outcome on the left-hand
#'   side, e.g. `Surv(time, status) ~ age + karno`.
#' @param data A data frame containing the formula variables.
#' @param model A learner id from `listlearners()$learner` (for example
#'   `"coxph"`, `"ranger"`, `"densemlp"`), or a custom learner: a list with
#'   `fit = function(formula, data, ...)` returning any fitted object and
#'   `predict = function(object, newdata, times)` returning an
#'   `nrow(newdata) x length(times)` matrix of survival probabilities. An
#'   optional `name` element labels it in printed output. A custom learner works
#'   with every verb.
#' @param spec A named list of engine arguments forwarded to the learner's
#'   fitting routine (for example `list(num.trees = 500)` for `"ranger"`).
#' @param seed Optional integer seed set immediately before fitting, for learners
#'   with a stochastic training step.
#' @param ... Not used; engine arguments passed here give an error. Use `spec`.
#'
#' @return An object of class `c("survalis_fit", "mlsurv_model")` with elements
#'   `model`, `learner`, `engine`, `formula`, `data`, `time`, `status`,
#'   `model_id`, and `spec`.
#'
#' @seealso [predict.survalis_fit()], [evaluate()], [tune()], [compare()],
#'   [interpret()], [estimate()], [listlearners()]
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran,
#'          model = "ranger", spec = list(num.trees = 200))
#' m
#'
#' # A custom learner: an exponential AFT model from survival::survreg
#' expo <- list(
#'   name = "exponential",
#'   fit = function(formula, data, ...) survival::survreg(formula, data, dist = "exponential"),
#'   predict = function(object, newdata, times) {
#'     rate <- 1 / predict(object, newdata, type = "response")
#'     exp(-outer(rate, times))
#'   }
#' )
#' m2 <- fit(Surv(time, status) ~ age + karno, veteran, model = expo)
#' predict(m2, veteran[1:3, ], times = c(90, 180))
#' @export
fit <- function(formula, data, model, spec = list(), seed = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  parsed <- tryCatch(.parse_surv_formula(formula, data), error = function(e) e)
  if (inherits(parsed, "error")) {
    stop("survalis::fit() could not read the formula: ", conditionMessage(parsed), "\n",
         "It needs a Surv(time, status) outcome on the left-hand side; ",
         "for non-survival outcomes use funcml::fit().", call. = FALSE)
  }

  lrn <- .resolve_learner(model)

  spec <- as.list(spec)
  dots <- list(...)
  if (length(dots)) {
    stop("fit(): engine arguments go in `spec = list(...)`; got unused argument(s) ",
         paste(sQuote(names(dots)), collapse = ", "), ".", call. = FALSE)
  }

  # Rows with a missing outcome or predictor are dropped once, here, so the
  # stored data, the printed counts and every later verb see the rows the
  # model was fitted on.
  data <- .complete_cases_df(data, all.vars(formula))
  if (!nrow(data)) stop("No complete rows for the variables in `formula`.", call. = FALSE)

  if (!is.null(seed)) set.seed(as.integer(seed))

  mlobj <- do.call(lrn$fit, c(list(formula = formula, data = data), spec))

  out <- .as_survalis_fit(mlobj, model_id = lrn$id, spec = spec,
                          formula = formula, data = data, parsed = parsed)
  if (lrn$custom) out$learner_def <- model
  out
}

#' @keywords internal
.null_default <- function(x, y) if (is.null(x)) y else x
`%||%` <- .null_default

#' Wrap a legacy mlsurv_model as a survalis_fit
#' @keywords internal
.as_survalis_fit <- function(mlobj, model_id, spec, formula, data, parsed) {
  if (!is.list(mlobj)) stop("fit_", model_id, "() did not return a model object.", call. = FALSE)
  mlobj$learner  <- .null_default(mlobj$learner, model_id)
  mlobj$formula  <- formula
  mlobj$data     <- data
  mlobj$time     <- parsed$time_col
  mlobj$status   <- parsed$status_col
  mlobj$model_id <- model_id
  mlobj$spec     <- spec
  mlobj$engine   <- .null_default(attr(mlobj, "engine"), NA_character_)
  attr(mlobj, "surv_parsed") <- parsed
  class(mlobj) <- unique(c("survalis_fit", setdiff(class(mlobj), "survalis_fit")))
  if (!"mlsurv_model" %in% class(mlobj)) class(mlobj) <- c(class(mlobj), "mlsurv_model")
  mlobj
}

#' 0/1 event vector of a fit's training data
#' @keywords internal
.fit_status_vector <- function(fit) {
  parsed <- attr(fit, "surv_parsed")
  if (is.null(parsed)) parsed <- .parse_surv_formula(fit$formula, fit$data)
  .recode_status_vec(parsed, fit$data)
}

#' @export
print.survalis_fit <- function(x, ...) {
  n <- nrow(x$data)
  st <- .fit_status_vector(x)
  cat(sprintf("<survalis_fit> model = %s%s\n", x$model_id,
              if (!is.na(x$engine)) sprintf("  (engine: %s)", x$engine) else ""))
  cat("  formula: ", deparse(x$formula), "\n", sep = "")
  cat(sprintf("  data:    %d obs, %d events (%.0f%%)\n",
              n, sum(st), 100 * mean(st)))
  if (length(x$spec)) {
    cat("  spec:    ",
        paste(names(x$spec), vapply(x$spec, function(v) paste(deparse(v), collapse = ""), ""),
              sep = " = ", collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.survalis_fit <- function(object, ...) {
  NextMethod()
}

#' @export
plot.survalis_fit <- function(x, times = NULL, group = NULL, ...) {
  S <- predict(x, newdata = x$data, times = times, type = "survival")
  plot_survmat(as.data.frame(S), times = attr(S, "times"), group = group, ...)
}
