#' Fit a survival learner
#'
#' The entry point of the survalis verb interface. `fit()` selects a learner by
#' id, passes engine arguments through `spec`, and returns a `survalis_fit`
#' object that the other verbs ([predict()], [evaluate()], [tune()], [interpret()],
#' [estimate()]) consume.
#'
#' `survalis_fit` extends the internal `mlsurv_model` contract, so every legacy
#' `predict_*()` / `compute_*()` helper keeps working on it unchanged.
#'
#' @param formula A survival formula with a `Surv()` outcome on the left-hand
#'   side, e.g. `Surv(time, status) ~ age + karno`.
#' @param data A data frame containing the formula variables.
#' @param model A learner id from `list_survlearners()$learner` (for example
#'   `"coxph"`, `"ranger"`, `"survdnn"`).
#' @param spec A named list of engine arguments forwarded to the underlying
#'   `fit_<model>()` (for example `list(num.trees = 500)` for `"ranger"`).
#' @param seed Optional integer seed set immediately before fitting, for learners
#'   with a stochastic training step.
#' @param ... Deprecated: engine arguments passed here are merged into `spec`
#'   with a one-time warning. Use `spec`.
#'
#' @return An object of class `c("survalis_fit", "mlsurv_model")` with elements
#'   `model`, `learner`, `engine`, `formula`, `data`, `time`, `status`,
#'   `model_id`, and `spec`.
#'
#' @seealso [predict.survalis_fit()], [evaluate()], [tune()], [compare()],
#'   [interpret()], [estimate()], [list_survlearners()]
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran,
#'          model = "ranger", spec = list(num.trees = 200))
#' m
#' @export
fit <- function(formula, data, model, spec = list(), seed = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  parsed <- tryCatch(.parse_surv_formula(formula, data), error = function(e) NULL)
  if (is.null(parsed)) {
    stop("survalis::fit() requires a Surv() outcome on the left-hand side.\n",
         "For non-survival outcomes use funcml::fit().", call. = FALSE)
  }

  known <- list_survlearners()$learner
  if (length(model) != 1L || !is.character(model) || !model %in% known) {
    stop("Unknown model ", sQuote(model %||% "NULL"),
         ". See list_survlearners() for the ", length(known), " available ids.",
         call. = FALSE)
  }

  spec <- as.list(spec)
  dots <- list(...)
  if (length(dots)) {
    rlang::warn(
      "fit(): pass engine arguments through `spec = list(...)`, not `...`. Merging them for now.",
      .frequency = "once", .frequency_id = "survalis-fit-dots"
    )
    spec <- utils::modifyList(spec, dots)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  fit_fun <- get(paste0("fit_", model), envir = asNamespace("survalis"))
  mlobj <- do.call(fit_fun, c(list(formula = formula, data = data), spec))

  .as_survalis_fit(mlobj, model_id = model, spec = spec,
                   formula = formula, data = data, parsed = parsed)
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
  if (is.null(mlobj$data)) mlobj$data <- data
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
  s <- fit$data[[parsed$status_col]]
  if (isTRUE(parsed$recode_status)) as.integer(s == parsed$event_value) else as.integer(s != 0)
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
