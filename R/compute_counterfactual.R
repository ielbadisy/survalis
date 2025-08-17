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





