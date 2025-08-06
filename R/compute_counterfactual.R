compute_counterfactual <- function(model, newdata, times, target_time,
                                   features_to_change = NULL, grid.size = 100,
                                   max.change = NULL, cost_penalty = 0.01) {
  if (nrow(newdata) != 1) {
    stop("newdata must contain exactly one row.")
  }

  # Auto-predictor from model
  learner <- model$learner
  pred_fun <- get(paste0("predict_", learner))

  # Convert newdata columns to factors if needed
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

    # Generate candidate values
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
      # Handle categorical
      if (is.factor(model$data[[feature]])) {
        newdata_mod[[feature]] <- factor(candidate_values[i], levels = levels(model$data[[feature]]))
      } else {
        newdata_mod[[feature]] <- candidate_values[i]
      }

      pred_mod <- pred_fun(model, newdata_mod, times = times)
      if (is.null(dim(pred_mod))) pred_mod <- matrix(pred_mod, nrow = 1)

      surv_mod_value <- pred_mod[1, idx_time]
      surv_gains[i] <- surv_mod_value - surv_orig_value

      # Cost calculation
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



library(ranger)
library(survival)

data(veteran)
veteran$trt <- factor(veteran$trt)
veteran$celltype <- factor(veteran$celltype)

mod_ranger <- fit_ranger(Surv(time, status) ~ trt + age + karno, data = veteran)

obs <- veteran[1, , drop = FALSE]

cf_result <- compute_counterfactual(
  model = mod_ranger,
  newdata = obs,
  times = c(50, 100, 150),
  target_time = 100,
  features_to_change = c("trt", "karno"),
  grid.size = 50,
  cost_penalty = 0.001
)

print(cf_result)


