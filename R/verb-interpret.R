#' Interpret a fitted survival learner
#'
#' A single entry point to the package's model-agnostic explanation methods. The
#' available `method` ids are derived from `list_interpretability_methods()`
#' (plus `"ice"`, a variant of `"pdp"`, and `"permute"`, an alias of
#' `"varimp"`).
#'
#' @param fit A `survalis_fit` (from [fit()]), a `survalis_tune`, or a legacy
#'   `mlsurv_model`.
#' @param method One of `"pdp"`, `"ice"`, `"ale"`, `"shap"`, `"shap_matrix"`,
#'   `"varimp"` / `"permute"`, `"interaction"`, `"counterfactual"`,
#'   `"surrogate"`, `"tree_surrogate"`, `"calibration"`, `"survcluster"`.
#' @param newdata Instances to explain (SHAP, ALE, counterfactual, surrogate,
#'   survcluster). Defaults to `data`.
#' @param data Background data for the method. Defaults to the training data.
#' @param times Evaluation grid; `NULL` uses [default_times()]. Methods that need
#'   a single time (`calibration`) take `eval_time` in `...` or the sole `times`
#'   value.
#' @param ... Method-specific arguments, e.g. `feature` (pdp / ice / ale),
#'   `target_time` (counterfactual), `eval_time` (calibration), `K` (survcluster).
#'
#' @return A `survalis_interpret` object wrapping the underlying `compute_*`
#'   result, with `method`, `times`, and a `plot()` method that dispatches to the
#'   paired `plot_*` helper.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
#'          spec = list(num.trees = 100))
#' interpret(m, "pdp", feature = "karno", times = c(90, 180))
#' interpret(m, "varimp", times = c(90, 180))
#' @export
interpret <- function(fit, method, newdata = NULL, data = NULL, times = NULL, ...) {
  if (inherits(fit, "survalis_tune")) fit <- fit$fit
  if (!inherits(fit, "mlsurv_model")) {
    stop("`fit` must come from fit() / tune().", call. = FALSE)
  }
  disp <- .interpret_dispatch()
  method <- match.arg(method, names(disp))
  entry <- disp[[method]]

  if (is.null(data)) data <- fit$data
  if (is.null(newdata)) newdata <- data

  parsed <- attr(fit, "surv_parsed")
  if (is.null(parsed)) parsed <- .parse_surv_formula(fit$formula, fit$data)
  status <- fit$data[[parsed$status_col]]
  status <- if (isTRUE(parsed$recode_status)) as.integer(status == parsed$event_value) else as.integer(status != 0)
  times <- .resolve_times(times, fit$data[[parsed$time_col]], status)

  dots <- list(...)
  result <- entry$run(fit, newdata = newdata, data = data, times = times, dots = dots)

  structure(
    list(method = method, result = result, plot_fn = entry$plot,
         times = times, feature = dots$feature, model = fit$learner),
    class = "survalis_interpret"
  )
}

#' Adaptors from interpret() method ids to compute_*() signatures
#' @keywords internal
.interpret_dispatch <- function() {
  need <- function(dots, nm, method) {
    if (is.null(dots[[nm]])) {
      stop("interpret(method = \"", method, "\") requires `", nm, "`.", call. = FALSE)
    }
    dots[[nm]]
  }
  drop_known <- function(dots, nms) dots[setdiff(names(dots), nms)]

  list(
    pdp = list(plot = "plot_pdp", run = function(fit, newdata, data, times, dots) {
      feature <- need(dots, "feature", "pdp")
      do.call(compute_pdp, c(list(model = fit, data = data, feature = feature,
        times = times, method = "pdp+ice"), drop_known(dots, "feature")))
    }),
    ice = list(plot = "plot_pdp", run = function(fit, newdata, data, times, dots) {
      feature <- need(dots, "feature", "ice")
      do.call(compute_pdp, c(list(model = fit, data = data, feature = feature,
        times = times, method = "ice"), drop_known(dots, c("feature", "method"))))
    }),
    ale = list(plot = "plot_ale", run = function(fit, newdata, data, times, dots) {
      feature <- need(dots, "feature", "ale")
      do.call(compute_ale, c(list(model = fit, newdata = newdata, feature = feature,
        times = times), drop_known(dots, "feature")))
    }),
    shap = list(plot = "plot_shap", run = function(fit, newdata, data, times, dots) {
      do.call(compute_shap, c(list(model = fit, newdata = newdata[1, , drop = FALSE],
        baseline_data = data, times = times), dots))
    }),
    shap_matrix = list(plot = "plot_shap_beeswarm", run = function(fit, newdata, data, times, dots) {
      do.call(compute_shap_matrix, c(list(model = fit, newdata = newdata,
        baseline_data = data, times = times), dots))
    }),
    varimp = list(plot = "plot_varimp", run = function(fit, newdata, data, times, dots) {
      do.call(compute_varimp, c(list(model = fit, times = times), dots))
    }),
    permute = list(plot = "plot_varimp", run = function(fit, newdata, data, times, dots) {
      do.call(compute_varimp, c(list(model = fit, times = times), dots))
    }),
    interaction = list(plot = "plot_interactions", run = function(fit, newdata, data, times, dots) {
      do.call(compute_interactions, c(list(model = fit, data = data, times = times), dots))
    }),
    counterfactual = list(plot = "plot_counterfactual", run = function(fit, newdata, data, times, dots) {
      tt <- need(dots, "target_time", "counterfactual")
      do.call(compute_counterfactual, c(list(model = fit, newdata = newdata,
        times = times, target_time = tt), drop_known(dots, "target_time")))
    }),
    surrogate = list(plot = "plot_surrogate", run = function(fit, newdata, data, times, dots) {
      do.call(compute_surrogate, c(list(model = fit, newdata = newdata,
        baseline_data = data), dots))
    }),
    tree_surrogate = list(plot = "plot_tree_surrogate", run = function(fit, newdata, data, times, dots) {
      do.call(compute_tree_surrogate, c(list(model = fit, data = data, times = times), dots))
    }),
    calibration = list(plot = "plot_calibration", run = function(fit, newdata, data, times, dots) {
      et <- dots$eval_time
      if (is.null(et)) {
        if (length(times) != 1L) {
          stop("interpret(method = \"calibration\") needs a single time: pass ",
               "`eval_time` or a length-1 `times`.", call. = FALSE)
        }
        et <- times
      }
      do.call(compute_calibration, c(list(model = fit, data = data,
        time = fit$time, status = fit$status, eval_time = et),
        drop_known(dots, "eval_time")))
    }),
    survcluster = list(plot = "plot_survcluster", run = function(fit, newdata, data, times, dots) {
      do.call(compute_survcluster, c(list(model = fit, newdata = newdata, times = times), dots))
    })
  )
}

#' @export
print.survalis_interpret <- function(x, ...) {
  cat(sprintf("<survalis_interpret> method = %s  |  model = %s\n", x$method, x$model))
  if (!is.null(x$feature)) cat("  feature:", x$feature, "\n")
  cat(sprintf("  grid: %d times in [%s, %s]\n",
              length(x$times), format(min(x$times)), format(max(x$times))))
  df <- tryCatch(as.data.frame(x$result), error = function(e) NULL)
  if (is.data.frame(df)) {
    cat("  $result:\n")
    print(utils::head(df, 5))
  } else {
    cat("  $result: <", class(x$result)[1], ">\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.survalis_interpret <- function(x, ...) {
  if (is.data.frame(x$result)) return(x$result)
  as.data.frame(x$result, ...)
}

#' @export
plot.survalis_interpret <- function(x, ...) {
  plot_fun <- get(x$plot_fn, envir = asNamespace("survalis"))
  args <- list(x$result, ...)
  if (!is.null(x$feature) && "feature" %in% names(formals(plot_fun))) {
    args$feature <- x$feature
  }
  do.call(plot_fun, args)
}
