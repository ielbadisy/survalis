compute_pdp <- function(model, predict_function, data, feature, times,
                        method = "pdp+ice", grid.size = 20) {

  if (!feature %in% names(data)) stop("Feature not found in data.")

  is_categorical <- is.factor(data[[feature]]) || is.character(data[[feature]])

  if (is_categorical) {
    grid_values <- sort(unique(data[[feature]]))
  } else {
    feature_range <- range(data[[feature]], na.rm = TRUE)
    grid_values <- seq(feature_range[1], feature_range[2], length.out = grid.size)
  }

  n_obs <- nrow(data)
  n_grid <- length(grid_values)
  expanded_data <- data[rep(seq_len(n_obs), each = n_grid), , drop = FALSE]
  expanded_data[[feature]] <- rep(grid_values, times = n_obs)

  surv_preds <- predict_function(model, newdata = expanded_data, times = times)
  if (is.null(dim(surv_preds))) surv_preds <- matrix(surv_preds, nrow = 1)

  surv_df <- as.data.frame(surv_preds)
  colnames(surv_df) <- paste0("t=", times)
  surv_df[[feature]] <- rep(grid_values, times = n_obs)
  surv_df$.id <- rep(seq_len(n_obs), each = n_grid)

  results_list <- list()

  if (method %in% c("ice", "pdp+ice")) {
    ice_long <- reshape(surv_df,
                        varying = paste0("t=", times),
                        v.names = "surv_prob",
                        timevar = "time",
                        times = times,
                        idvar = c(".id", feature),
                        direction = "long")
    ice_long$type <- "ice"
    rownames(ice_long) <- NULL
    results_list$ice <- ice_long
  }

  if (method %in% c("pdp", "pdp+ice")) {
    # Compute mean survival for each time and feature value
    means <- aggregate(surv_df[, paste0("t=", times)],
                       by = list(surv_df[[feature]]),
                       FUN = mean)
    colnames(means)[1] <- feature
    pdp_long <- reshape(means,
                        varying = paste0("t=", times),
                        v.names = "surv_prob",
                        timevar = "time",
                        times = times,
                        idvar = feature,
                        direction = "long")
    pdp_long$type <- "pdp"
    pdp_long$.id <- NA
    rownames(pdp_long) <- NULL
    results_list$pdp <- pdp_long
  }

  final_result <- do.call(rbind, results_list)
  final_result$time <- as.numeric(as.character(final_result$time))

  pdp_integrated <- NULL
  if ("pdp" %in% final_result$type) {
    integrate_fun <- function(x, y) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      sum((head(y, -1) + tail(y, -1)) / 2 * diff(x))
    }

    split_pdp <- split(final_result[final_result$type == "pdp", ], final_result[final_result$type == "pdp", feature])
    pdp_integrated <- do.call(rbind, lapply(split_pdp, function(df) {
      int_surv <- integrate_fun(df$time, df$surv_prob)
      f_val <- unique(df[[feature]])
      out <- data.frame(integrated_surv = int_surv, stringsAsFactors = FALSE)
      out[[feature]] <- f_val
      out
    }))
    rownames(pdp_integrated) <- NULL
    pdp_integrated <- pdp_integrated[, c(feature, "integrated_surv")]
  }

  list(
    results = final_result,
    pdp_integrated = pdp_integrated
  )
}


mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)

pdp_ice_result <- compute_pdp(
  model = mod_ranger,
  predict_function = predict_ranger,
  data = veteran,
  feature = "celltype",  # try "celltype" too
  times = c(100, 200, 300),
  method = "pdp+ice"
)

