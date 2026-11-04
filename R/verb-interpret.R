#' Interpret a fitted survival learner
#'
#' A single entry point to the package's model-agnostic explanation and
#' diagnostic methods. `list_interpretability_methods()` lists the `method` ids
#' with a one-line description of each.
#'
#' @param fit A `survalis_fit` (from [fit()]), a `survalis_tune`, or a legacy
#'   `mlsurv_model`.
#' @param method One of `"pdp"`, `"ice"`, `"ale"`, `"shap"`, `"shap_matrix"`,
#'   `"varimp"` / `"permute"`, `"interaction"`, `"counterfactual"`,
#'   `"surrogate"`, `"tree_surrogate"`, `"calibration"`, `"roc"`, `"dca"`,
#'   `"survcluster"`.
#' @param newdata Instances to explain (SHAP, ALE, counterfactual, surrogate,
#'   survcluster) or to assess (`roc`, `dca`; these need the outcome columns).
#'   Defaults to `data`, so `roc` and `dca` are apparent unless you pass
#'   held-out data.
#' @param data Background data for the method. Defaults to the training data.
#' @param times Evaluation grid; `NULL` picks one from the data quantiles.
#'   Methods that need a single time (`calibration`, `roc`, `dca`) take
#'   `eval_time` in `...` or the sole `times` value.
#' @param ... Method-specific arguments, e.g. `feature` (pdp / ice / ale),
#'   `target_time` (counterfactual, surrogate), `eval_time` (calibration, roc,
#'   dca), `thresholds` (dca), `K` (survcluster). For `"shap"`, the aggregation
#'   over time (`"meanabs"` or `"integral"`) is passed as `shap_method`.
#'
#' @return A `survalis_interpret` object holding the method's `result`, with
#'   `method`, `times`, and a `plot()` method. `as.data.frame()` returns the
#'   result frame. For `"survcluster"` it is in long format, one row per
#'   observation and time: `id` (row of `newdata`), `cluster`, `medoid` (`TRUE`
#'   for the cluster's medoid curve), `time` and predicted survival `surv`.
#'   For `"roc"` and `"dca"` it is the curve with an `eval_time` column (and
#'   `auc` for `"roc"`); for `"tree_surrogate"` one row per time and feature
#'   with the number of `splits` on it and the tree's `r_squared`. For
#'   `"pdp"` / `"ice"` it is one row per curve point (`feature`,
#'   `feature_value`, `id` (`NA` on PDP rows), `time`, `surv`, `type`), for
#'   `"ale"` one row per feature value and time (`feature`, `feature_value`,
#'   `time`, `ale`), and for `"calibration"` the binned table with its
#'   `eval_time`. For `"pdp"` and `"ale"`, `as.data.frame(x, what =
#'   "integrated")` returns the curve integrated over time instead (`feature`,
#'   `feature_value`, `integrated`).
#'
#' @seealso For scalar summaries use [evaluate()]: `metrics = "ece"` summarises
#'   calibration and `metrics = "auc"` summarises the ROC curve.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
#'          spec = list(num.trees = 100))
#' interpret(m, "pdp", feature = "karno", times = c(90, 180))
#' interpret(m, "varimp", times = c(90, 180))
#' plot(interpret(m, "roc", eval_time = 180))
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
  one_time <- function(dots, times, method) {
    et <- dots$eval_time
    if (is.null(et)) {
      if (length(times) != 1L) {
        stop("interpret(method = \"", method, "\") needs a single time: pass ",
             "`eval_time` or a length-1 `times`.", call. = FALSE)
      }
      et <- times
    }
    et
  }
  outcome_and_pred <- function(fit, newdata, et) {
    parsed <- attr(fit, "surv_parsed")
    if (is.null(parsed)) parsed <- .parse_surv_formula(fit$formula, fit$data)
    need <- c(parsed$time_col, parsed$status_col)
    if (!all(need %in% names(newdata))) {
      stop("`newdata` must contain the outcome columns ",
           paste(sQuote(need), collapse = ", "), ".", call. = FALSE)
    }
    newdata <- .complete_cases_df(newdata, all.vars(fit$formula))
    y <- survival::Surv(newdata[[parsed$time_col]], .recode_status_vec(parsed, newdata))
    p <- .finalize_survmat(.predict_fun(fit)(fit, newdata = newdata, times = et), times = et)
    list(y = y, p = p)
  }

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
      # `method` names interpret()'s own dispatch argument, so compute_shap()'s
      # aggregation method (also called `method`) is reachable as `shap_method`.
      if (!is.null(dots$shap_method)) {
        dots$method <- dots$shap_method
        dots$shap_method <- NULL
      }
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
      tt <- need(dots, "target_time", "surrogate")
      do.call(compute_surrogate, c(list(model = fit, newdata = newdata,
        baseline_data = data, times = times, target_time = tt),
        drop_known(dots, "target_time")))
    }),
    tree_surrogate = list(plot = "plot_tree_surrogate", run = function(fit, newdata, data, times, dots) {
      do.call(compute_tree_surrogate, c(list(model = fit, data = data, times = times), dots))
    }),
    calibration = list(plot = "plot_calibration", run = function(fit, newdata, data, times, dots) {
      et <- one_time(dots, times, "calibration")
      do.call(compute_calibration, c(list(model = fit, data = data,
        time = fit$time, status = fit$status, eval_time = et),
        drop_known(dots, "eval_time")))
    }),
    roc = list(plot = "plot_roc", run = function(fit, newdata, data, times, dots) {
      et <- one_time(dots, times, "roc")
      op <- outcome_and_pred(fit, newdata, et)
      roc_survmat(op$y, predicted = op$p, t_star = et)
    }),
    dca = list(plot = "plot_dca", run = function(fit, newdata, data, times, dots) {
      et <- one_time(dots, times, "dca")
      op <- outcome_and_pred(fit, newdata, et)
      args <- list(op$y, predicted = op$p, t_star = et)
      if (!is.null(dots$thresholds)) args$thresholds <- dots$thresholds
      do.call(dca_survmat, args)
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
  cat("  $result:\n")
  print(utils::head(as.data.frame(x), 5))
  invisible(x)
}

#' @export
as.data.frame.survalis_interpret <- function(x, row.names = NULL, optional = FALSE,
                                            what = c("curves", "integrated"), ...) {
  what <- match.arg(what)
  res <- x$result
  if (what == "integrated") {
    if (!x$method %in% c("pdp", "ale")) {
      stop("`what = \"integrated\"` applies to interpret(method = \"pdp\") and ",
           "\"ale\" only.", call. = FALSE)
    }
    return(.integrated_frame(res, x$method, x$feature))
  }
  if (is.data.frame(res)) {
    # plain data.frame (varimp comes back as a basetable); attributes are kept
    class(res) <- "data.frame"
    return(res)
  }
  switch(x$method,
    pdp = ,
    ice = .pdp_long(res),
    ale = .ale_long(res, x$feature),
    calibration = data.frame(eval_time = res$eval_time, res$calibration_table),
    survcluster = .survcluster_long(res),
    roc = data.frame(eval_time = res$t_star, res$curve, auc = res$auc),
    dca = data.frame(eval_time = res$t_star, res$curve),
    tree_surrogate = .tree_surrogate_long(res),
    stop("No data frame form for interpret(method = \"", x$method, "\").",
         call. = FALSE)
  )
}

#' Long format of a pdp / ice result: one row per curve point
#' @noRd
.pdp_long <- function(res) {
  r <- res$results
  feature <- names(r)[1]
  out <- data.frame(feature = feature, feature_value = r[[1]], id = r$.id,
                    time = r$time, surv = r$surv_prob, type = r$type,
                    stringsAsFactors = FALSE)
  out
}

#' Integrated (time-averaged) curve of a pdp or ale result
#' @noRd
.integrated_frame <- function(res, method, feature) {
  integ <- if (method == "pdp") res$pdp_integrated else res$integrated
  if (!is.data.frame(integ) || !nrow(integ)) {
    stop("This result has no integrated curve (it needs a PDP over more than ",
         "one time).", call. = FALSE)
  }
  feat <- if (method == "pdp") names(integ)[1] else feature
  data.frame(feature = if (is.null(feat)) NA_character_ else feat,
             feature_value = integ[[1]],
             integrated = integ[[2]], stringsAsFactors = FALSE)
}

#' Long format of an ale result: one row per feature value and time
#' @noRd
.ale_long <- function(res, feature) {
  ale <- res$ale
  tcols <- setdiff(names(ale), "feature_value")
  times <- as.numeric(sub("^t=", "", tcols))
  nv <- nrow(ale)
  out <- data.frame(
    feature = if (is.null(feature)) NA_character_ else feature,
    feature_value = rep(ale$feature_value, times = length(tcols)),
    time = rep(times, each = nv),
    ale = unlist(ale[tcols], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  out
}

#' Long format of a survcluster result: one row per observation and time
#' @noRd
.survcluster_long <- function(res) {
  S <- res$survmat
  n <- nrow(S)
  nt <- length(res$times)
  cl <- res$fit$clusters
  med <- seq_len(n) %in% res$fit$medoid_indices
  data.frame(
    id = rep(seq_len(n), times = nt),
    cluster = rep(cl, times = nt),
    medoid = rep(med, times = nt),
    time = rep(res$times, each = n),
    surv = as.vector(S),
    stringsAsFactors = FALSE
  )
}

#' Long format of a tree_surrogate result: one row per time and feature
#' @noRd
.tree_surrogate_long <- function(res) {
  parts <- lapply(res$results, function(r) {
    sc <- r$split_count
    if (!length(sc)) {
      return(data.frame(time = r$time, feature = NA_character_, splits = 0L,
                        r_squared = r$r_squared, stringsAsFactors = FALSE))
    }
    data.frame(time = r$time, feature = names(sc), splits = as.integer(sc),
               r_squared = r$r_squared, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
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
