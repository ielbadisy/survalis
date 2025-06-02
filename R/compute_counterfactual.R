counterfactual_surv_v3 <- function(model, predict_function, newdata, times, target_time,
                                   features_to_change = NULL, grid.size = 100,
                                   max.change = NULL, cost_penalty = 0.01) {
  requireNamespace("data.table")

  if (nrow(newdata) != 1) {
    stop("newdata must contain exactly one row.")
  }

  feature_names <- setdiff(colnames(newdata), c(model$time, model$status))

  if (!is.null(features_to_change)) {
    feature_names <- intersect(feature_names, features_to_change)
  }

  if (length(feature_names) == 0) {
    stop("No valid features to change.")
  }

  idx_time <- which.min(abs(times - target_time))

  surv_orig <- predict_function(model, newdata, times = times)
  if (is.null(dim(surv_orig))) surv_orig <- matrix(surv_orig, nrow = 1)
  surv_orig_value <- surv_orig[1, idx_time]

  results <- list()

  for (feature in feature_names) {
    value_orig <- newdata[[feature]]

    # grid generation
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
      candidate_values <- unique(newdata[[feature]])
      candidate_values <- setdiff(candidate_values, value_orig)
      if (length(candidate_values) == 0) next
    } else {
      next
    }

    n_candidates <- length(candidate_values)
    surv_gains <- numeric(n_candidates)
    change_costs <- numeric(n_candidates)

    for (i in seq_along(candidate_values)) {
      newdata_mod <- newdata
      newdata_mod[[feature]] <- candidate_values[i]

      pred_mod <- predict_function(model, newdata_mod, times = times)
      if (is.null(dim(pred_mod))) pred_mod <- matrix(pred_mod, nrow = 1)

      surv_mod_value <- pred_mod[1, idx_time]
      surv_gains[i] <- surv_mod_value - surv_orig_value

      # Cost
      if (is.numeric(value_orig)) {
        change_costs[i] <- abs(as.numeric(candidate_values[i]) - as.numeric(value_orig))
      } else {
        # categorical change cost = 1 (switch)
        change_costs[i] <- 1
      }
    }

    # Penalized objective
    penalized_gain <- surv_gains - cost_penalty * change_costs

    best_idx <- which.max(penalized_gain)

    results[[feature]] <- data.frame(
      feature = feature,
      original_value = value_orig,
      suggested_value = candidate_values[best_idx],
      survival_gain = surv_gains[best_idx],
      change_cost = change_costs[best_idx],
      penalized_gain = penalized_gain[best_idx],
      stringsAsFactors = FALSE
    )
  }

  final_result <- do.call(rbind, results)
  rownames(final_result) <- NULL
  final_result
}



#***********

predict_ranger <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "ranger")
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  pred <- predict(object$model, data = newdata, type = "response")

  surv_mat <- pred$survival
  model_times <- pred$unique.death.times

  if (is.null(times)) {
    times <- model_times
  }

  # Ensure matrix shape even if only one row
  if (is.null(dim(surv_mat))) {
    surv_mat <- matrix(surv_mat, nrow = 1)
  }

  # Interpolate survival probabilities
  surv_probs <- apply(surv_mat, 1, function(row_surv) {
    stats::approx(
      x = model_times,
      y = row_surv,
      xout = times,
      method = "linear",
      rule = 2
    )$y
  })

  # Ensure correct orientation
  if (length(times) == 1) {
    surv_probs <- matrix(surv_probs, ncol = 1)
  } else {
    surv_probs <- t(surv_probs)
  }

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  as.data.frame(surv_probs)
}


#test

cf_result <- counterfactual_surv_v3(
  model = mod_ranger,
  predict_function = predict_ranger,
  newdata = sim_data[1, , drop = FALSE],
  times = c(50, 100, 150),
  target_time = 300,
  features_to_change = c("trt", "karno"),
  grid.size = 100,
  max.change = list(age = 20, karno = 30),
  cost_penalty = 0.001
)

print(cf_result)


