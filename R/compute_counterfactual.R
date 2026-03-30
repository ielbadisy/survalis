#' Compute individual counterfactual changes to increase survival
#'
#' For a single individual (one-row \code{newdata}), propose feature changes
#' that maximize survival probability at a target time, subject to optional
#' per-feature change costs and bounds inferred from the training data in
#' \code{model$data}.
#'
#' @param model A fitted survival model (e.g., an \code{mlsurv_model}) that was
#'   trained with a data frame stored in \code{model$data}, and (ideally) names
#'   of the outcome variables in \code{model$time} and \code{model$status}.
#' @param newdata A data frame with exactly one row representing the individual.
#' @param times Numeric vector of time points used for prediction. Required unless
#'   the model's predict function can infer times; used together with \code{target_time}.
#' @param target_time Numeric scalar time at which to optimize survival. If missing
#'   and \code{times} is provided, the median of \code{times} is used. If \code{times}
#'   is missing but \code{target_time} is provided, \code{times} is set to \code{target_time}.
#' @param features_to_change Optional character vector of feature names allowed to change.
#'   Defaults to all predictors (non-outcome columns).
#' @param grid.size Integer grid size for numeric features (default 100).
#' @param max.change Optional named list of numeric bounds for per-feature absolute change,
#'   e.g., \code{list(age = 5)}.
#' @param cost_penalty Numeric penalty weight applied to magnitude of change (default 0.01).
#'
#' @return A \code{data.frame} with one row per feature considered and columns:
#'   \code{feature}, \code{original_value}, \code{suggested_value},
#'   \code{survival_gain}, \code{change_cost}, \code{penalized_gain}.
#'
#' @details
#' For each candidate feature, the function sweeps over plausible values
#' (numeric grid between observed min/max; all other levels for categorical),
#' predicts survival at \code{target_time}, and reports the best penalized gain
#' relative to the original value. Survival predictions are obtained via a
#' corresponding \code{predict_*} function inferred from \code{model$learner}
#' (e.g., \code{predict_coxph} for \code{learner = "coxph"}).
#'
#' @examples
#' # Uses only {survival} (already in Imports) and your internal `veteran` dataset
#' veteran$A <- veteran$trt
#' mod <- fit_coxph(survival::Surv(time, status) ~ A + age + karno, data = veteran)
#'
#' cf <- compute_counterfactual(
#'   model = mod,
#'   newdata = veteran[1, , drop = FALSE],   # exactly one individual
#'   times = c(50, 100, 150),
#'   target_time = 100,
#'   features_to_change = c("A", "age", "karno"),
#'   cost_penalty = 0.01
#' )
#' head(cf)
#'
#' @export

compute_counterfactual <- function(model, newdata, times, target_time,
                                   features_to_change = NULL, grid.size = 100,
                                   max.change = NULL, cost_penalty = 0.01) {
  if (nrow(newdata) != 1) {
    stop("newdata must contain exactly one row.")
  }

  # infer the predict function 
  learner <- model$learner
  pred_fun <- get(paste0("predict_", learner))

  # convert newdata columns to factors if needed
  for (col in names(newdata)) {
    if (is.factor(model$data[[col]]) && !is.factor(newdata[[col]])) {
      newdata[[col]] <- factor(newdata[[col]], levels = levels(model$data[[col]]))
    }
  }

  feature_names <- setdiff(colnames(newdata), c(model$time, model$status))

  if (!is.null(features_to_change)) {
    feature_names <- intersect(feature_names, features_to_change)
  }

  if (length(feature_names) == 0) {
    stop("No valid features to change.")
  }

  idx_time <- which.min(abs(times - target_time))

  surv_orig <- pred_fun(model, newdata, times = times)
  if (is.null(dim(surv_orig))) surv_orig <- matrix(surv_orig, nrow = 1)
  surv_orig_value <- surv_orig[1, idx_time]

  results <- list()

  for (feature in feature_names) {
    value_orig <- newdata[[feature]]

    # generate candidate values
    if (is.numeric(value_orig)) {
      observed_min <- min(model$data[[feature]], na.rm = TRUE)
      observed_max <- max(model$data[[feature]], na.rm = TRUE)
      lower <- observed_min
      upper <- observed_max
      if (!is.null(max.change) && feature %in% names(max.change)) {
        lower <- max(value_orig - max.change[[feature]], lower)
        upper <- min(value_orig + max.change[[feature]], upper)
      }
      candidate_values <- seq(lower, upper, length.out = grid.size)
    } else if (is.factor(value_orig) || is.character(value_orig)) {
      levels_full <- if (is.factor(model$data[[feature]])) {
        levels(model$data[[feature]])
      } else {
        unique(model$data[[feature]])
      }
      candidate_values <- setdiff(levels_full, as.character(value_orig))
      if (length(candidate_values) == 0) next
    } else {
      next
    }

    surv_gains <- numeric(length(candidate_values))
    change_costs <- numeric(length(candidate_values))

    for (i in seq_along(candidate_values)) {
      newdata_mod <- newdata
      # handle categorical
      if (is.factor(model$data[[feature]])) {
        newdata_mod[[feature]] <- factor(candidate_values[i], levels = levels(model$data[[feature]]))
      } else {
        newdata_mod[[feature]] <- candidate_values[i]
      }

      pred_mod <- pred_fun(model, newdata_mod, times = times)
      if (is.null(dim(pred_mod))) pred_mod <- matrix(pred_mod, nrow = 1)

      surv_mod_value <- pred_mod[1, idx_time]
      surv_gains[i] <- surv_mod_value - surv_orig_value

      # cost calculation
      if (is.numeric(value_orig)) {
        change_costs[i] <- abs(as.numeric(candidate_values[i]) - as.numeric(value_orig))
      } else {
        change_costs[i] <- 1  # flat cost for categorical
      }
    }

    penalized_gain <- surv_gains - cost_penalty * change_costs
    best_idx <- which.max(penalized_gain)

    results[[feature]] <- data.frame(
      feature = feature,
      original_value = as.character(value_orig),
      suggested_value = if (is.numeric(value_orig)) {
        as.character(round(candidate_values[best_idx], 4))
      } else {
        as.character(candidate_values[best_idx])
      },
      survival_gain = round(surv_gains[best_idx], 4),
      change_cost = round(change_costs[best_idx], 4),
      penalized_gain = round(penalized_gain[best_idx], 4),
      stringsAsFactors = FALSE
    )
  }

  final_result <- do.call(rbind, results)
  rownames(final_result) <- NULL
  return(final_result)
}

#' Plot Counterfactual Recommendations
#'
#' Visualizes counterfactual feature changes returned by
#' \code{\link{compute_counterfactual}}, ranking recommendations by either raw
#' survival gain or penalized gain.
#'
#' @param counterfactual_df A data frame returned by \code{compute_counterfactual()}.
#' @param metric Character scalar; one of \code{"penalized_gain"} (default) or
#'   \code{"survival_gain"} indicating which column to plot on the x-axis.
#' @param top_n Optional integer limiting the plot to the top \code{n} features
#'   after sorting by \code{metric}.
#' @param include_negative Logical; if \code{FALSE} (default), only positive-gain
#'   recommendations are shown.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' # cf <- compute_counterfactual(
#' #   mod, newdata = veteran[1, , drop = FALSE], times = 100, target_time = 100
#' # )
#' # plot_counterfactual(cf)
#' @export
plot_counterfactual <- function(counterfactual_df,
                                metric = c("penalized_gain", "survival_gain"),
                                top_n = NULL,
                                include_negative = FALSE) {
  metric <- match.arg(metric)

  required_cols <- c(
    "feature", "original_value", "suggested_value",
    "survival_gain", "change_cost", "penalized_gain"
  )
  missing_cols <- setdiff(required_cols, names(counterfactual_df))
  if (length(missing_cols)) {
    stop(
      "`counterfactual_df` must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  plot_df <- as.data.frame(counterfactual_df, stringsAsFactors = FALSE)
  plot_df[[metric]] <- as.numeric(plot_df[[metric]])

  if (!isTRUE(include_negative)) {
    positive_df <- plot_df[plot_df[[metric]] > 0, , drop = FALSE]
    if (nrow(positive_df)) {
      plot_df <- positive_df
    }
  }

  plot_df <- plot_df[order(plot_df[[metric]], decreasing = TRUE), , drop = FALSE]
  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (!is.finite(top_n) || top_n < 1) {
      stop("`top_n` must be a positive integer.")
    }
    plot_df <- utils::head(plot_df, top_n)
  }

  plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
  plot_df$change_label <- paste0(plot_df$original_value, " -> ", plot_df$suggested_value)
  plot_df$direction <- ifelse(plot_df[[metric]] >= 0, "increase", "decrease")
  plot_df$hjust <- ifelse(plot_df[[metric]] >= 0, -0.05, 1.05)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[metric]], y = feature, fill = direction)
  ) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = change_label, hjust = hjust),
      size = 3.5
    ) +
    ggplot2::scale_fill_manual(values = c(increase = "steelblue", decrease = "tomato")) +
    ggplot2::labs(
      title = "Counterfactual feature recommendations",
      x = if (metric == "penalized_gain") "Penalized survival gain" else "Survival gain",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::xlim(
      min(0, min(plot_df[[metric]], na.rm = TRUE)),
      max(plot_df[[metric]], na.rm = TRUE) * 1.25
    )
}
