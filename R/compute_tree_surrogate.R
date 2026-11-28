
#' Compute Tree-Based Surrogate Model for Survival Predictions
#'
#' @description
#' Fits decision tree surrogate models to approximate the predictions of a fitted survival model
#' at one or more evaluation times. This allows users to gain interpretable, rule-based
#' approximations of complex survival models.
#'
#' @param model A fitted survival model object created with a `fit_*()` function.
#'              Must include a valid `learner` field corresponding to a `predict_*()` method.
#' @param data A data frame containing predictor variables (and optional survival outcome columns).
#' @param times A numeric vector of evaluation times at which to approximate model predictions.
#'              Must contain at least one value.
#' @param minsplit Minimum number of observations required to attempt a split in the surrogate tree.
#'                 Passed to `rpart::rpart.control()`.
#' @param cp Complexity parameter for the surrogate tree. Passed to `rpart::rpart.control()`.
#' @param newdata Optional data frame of other rows on which to measure the tree's fidelity
#'   out of sample (`r_squared_test`).
#' @param holdout Fraction of the rows of `data` held out at random (before fitting the tree)
#'   to measure fidelity out of sample when `newdata` is not given; `0` (default) holds none.
#'
#' @details
#' For each evaluation time, the function:
#' 1. Predicts survival probabilities from the fitted model.
#' 2. Excludes survival outcome columns (`time`, `status`, `event`) from the predictors.
#' 3. Fits a decision tree to approximate the predicted probabilities.
#' 4. Computes the fidelity R^2 = 1 - SSE/SST of the tree's predictions against the model's,
#'    on the rows the tree was grown on (`r_squared`) and, with `newdata` or `holdout`, on
#'    other rows (`r_squared_test`, SST about their own mean).
#' 5. Counts the number of splits per feature.
#'
#' If multiple `times` are provided, results are stored for each time point.
#'
#' @return
#' An object of class `"tree_surrogate"`, containing:
#' - `times`: the evaluation times.
#' - `results`: a list with one element per time, each containing:
#'   - `tree`: the fitted `rpart` object.
#'   - `r_squared`: the in-sample fidelity R^2 of the tree vs. the model's predictions.
#'   - `r_squared_test`: the out-of-sample fidelity R^2 (`NA` without `newdata` or `holdout`).
#'   - `split_count`: a table of feature split counts.
#' - `dynamic`: logical indicating if more than one time was used.
#'
#' @examples
#' \dontrun{
#' mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' tree_ranger <- compute_tree_surrogate(
#'   model = mod_ranger,
#'   data = veteran,
#'   times = c(100, 200, 300)
#' )
#' }
#'
#' @noRd

compute_tree_surrogate <- function(model, data, times, minsplit = 10, cp = 0.01, newdata = NULL,
                                   holdout = 0) {
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Please install 'rpart'.")
  if (!requireNamespace("partykit", quietly = TRUE)) stop("Please install 'partykit'.")

  if (length(times) < 1) stop("At least one evaluation time must be specified.")

  # infer prediction function from model
  if (is.null(model$learner) && !is.function(model$predict_fun)) {
    stop("`model` must come from fit().", call. = FALSE)
  }
  predict_function <- .predict_fun(model)

  is_dynamic <- length(times) > 1
  if (is.null(newdata) && holdout > 0) {
    i <- sample.int(nrow(data), round(holdout * nrow(data)))
    newdata <- data[i, , drop = FALSE]
    data <- data[-i, , drop = FALSE]
  }
  covariates <- .model_features(model, names(data))
  r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
  surv_at <- function(d, time) {
    preds <- predict_function(model, newdata = d, times = time)
    col_match <- grep(paste0("^t=*", time, "$"), colnames(preds))
    if (length(col_match) == 1) return(preds[[col_match]])
    warning(sprintf("Could not find a single matching column for time = %s. Using first column.", time))
    if (is.data.frame(preds)) preds[[1]] else preds[, 1]
  }

  results <- lapply(times, function(time) {
    y <- surv_at(data, time)
    dat <- cbind(data[, covariates, drop = FALSE], surv_prob = y)

    # Fit surrogate tree
    tree_fit <- rpart::rpart(
      surv_prob ~ .,
      data = dat,
      control = rpart::rpart.control(minsplit = minsplit, cp = cp)
    )

    r2_test <- if (is.null(newdata)) NA_real_ else
      r2(surv_at(newdata, time), predict(tree_fit, newdata = newdata[, covariates, drop = FALSE]))

    split_count <- table(tree_fit$frame$var[tree_fit$frame$var != "<leaf>"])

    list(
      time = time,
      tree = tree_fit,
      r_squared = r2(y, predict(tree_fit, newdata = dat)),
      r_squared_test = r2_test,
      split_count = split_count
    )
  })

  names(results) <- paste0("time_", times)

  structure(list(
    times = times,
    results = results,
    dynamic = is_dynamic
  ), class = "tree_surrogate")
}


#' Plot Tree-Based Surrogate Models or Feature Importances
#'
#' @description
#' Visualizes the results of a `tree_surrogate` object returned by
#' `compute_tree_surrogate()`. Can display either the fitted surrogate trees
#' or aggregated feature importance based on split counts.
#'
#' @param tree_surrogate An object of class `"tree_surrogate"`.
#' @param type Character string indicating the type of plot:
#'   - `"tree"`: displays the surrogate decision tree(s) for each time.
#'   - `"importance"`: displays a bar plot of top features ranked by split count.
#' @param top_n Integer, the number of top features to display in the importance plot.
#'              Ignored if `type = "tree"`.
#' @param title Plot title for `type = "importance"`. If missing, an
#'   automatically generated title is used. Pass `NULL` to omit the title
#'   entirely (e.g., for journals requiring caption-only figures). Ignored
#'   for `type = "tree"`, which titles each per-time tree panel individually.
#'
#' @details
#' - If `type = "tree"`, a separate tree diagram is produced for each evaluation time.
#' - If `type = "importance"`, feature split counts are summed across all times
#'   and plotted as a bar chart.
#'
#' Requires the `partykit` package for tree plotting and `ggplot2` for importance plotting.
#'
#' @return
#' A plot object (for `"importance"`) or printed tree diagrams (for `"tree"`).
#'
#' @examples
#' \dontrun{
#' plot_tree_surrogate(tree_ranger, type = "tree")
#' plot_tree_surrogate(tree_ranger, type = "importance", top_n = 5)
#' }
#'
#' @noRd

plot_tree_surrogate <- function(tree_surrogate, type = c("tree", "importance"), top_n = 10, title) {
  if (!requireNamespace("partykit", quietly = TRUE)) stop("Please install 'partykit'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  type <- match.arg(type)

  results <- tree_surrogate$results

  if (type == "tree") {
    for (res in results) {
      party_tree <- partykit::as.party(res$tree)
      plot(
        party_tree,
        main = paste0("Surrogate Tree at time = ", res$time, ", R2 = ", round(res$r_squared, 2),
                      if (!is.na(.null_default(res$r_squared_test, NA)))
                        paste0(" (out of sample ", round(res$r_squared_test, 2), ")"))
      )
    }
  } else if (type == "importance") {
    importances <- lapply(results, function(res) {
      sc <- res$split_count
      data.frame(feature = names(sc), count = as.integer(sc), time = res$time)
    })
    all_importance <- basetable::rbindfill(importances, fill = TRUE)
    avg_importance <- basetable::aggregate(all_importance, by = "feature",
                                           value = "count", fun = sum)
    avg_importance <- basetable::orderrows(avg_importance, by = "count", decreasing = TRUE)
    avg_importance <- avg_importance[seq_len(min(top_n, nrow(avg_importance))), , drop = FALSE]

    if (missing(title)) title <- "Top Features by Split Count (All Times)"
    p <- ggplot2::ggplot(avg_importance, ggplot2::aes(x = reorder(feature, count), y = count)) +
      ggplot2::geom_bar(stat = "identity", fill = .survalis_palette[1]) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Feature",
        y = "Split Count"
      ) +
      theme_survalis()
    if (!is.null(title)) p <- p + ggplot2::labs(title = title)
    p
  }
}




