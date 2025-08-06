compute_pdp <- function(model, data, feature, times,
                        method = "pdp+ice", grid.size = 20) {
  # Automatically determine prediction function
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model has a valid `learner` and was created using `fit_*()`.")
  }
  predict_function <- get(paste0("predict_", model$learner))

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
  if (is.null(dim(surv_preds))) surv_preds <- matrix(surv_preds, ncol = 1)
  if (is.null(colnames(surv_preds))) colnames(surv_preds) <- paste0("t=", times)

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
    means <- aggregate(surv_df[, paste0("t=", times), drop = FALSE],
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
  if ("pdp" %in% final_result$type && length(unique(final_result$time)) > 1) {
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




plot_pdp <- function(pdp_ice_output, feature,
                     method = "pdp+ice", ids = NULL,
                     which = c("per_time", "integrated"),
                     alpha_ice = 0.2, smooth = FALSE) {
  requireNamespace("ggplot2")
  requireNamespace("data.table")

  which <- match.arg(which)
  results <- data.table::as.data.table(pdp_ice_output$results)
  feature_data <- results[[feature]]
  is_categorical <- is.factor(feature_data) || is.character(feature_data)

  if (which == "integrated") {

    if (is.null(pdp_ice_output$pdp_integrated)) stop("No integrated PDP found.")

    if (!"integrated_surv" %in% names(pdp_ice_output$pdp_integrated)) stop("Invalid integrated PDP.")

    ice_data <- results[type == "ice"]
    ice_integrated <- ice_data[, .(
      integrated_surv = sum((head(surv_prob, -1) + tail(surv_prob, -1)) / 2 * diff(head(time, -1))),
      time_n = .N
    ), by = c(".id", feature)]

    pdp_integrated <- pdp_ice_output$pdp_integrated

    if (is_categorical) {
      p <- ggplot2::ggplot() +
        ggplot2::geom_boxplot(data = ice_integrated, ggplot2::aes_string(x = feature, y = "integrated_surv"),
                              alpha = alpha_ice, fill = "pink") +
        ggplot2::geom_point(data = pdp_integrated, ggplot2::aes_string(x = feature, y = "integrated_surv"),
                            shape = 21, size = 3, fill = "black") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste("Integrated PDP with ICE Boxplot for", feature),
          x = feature,
          y = "Integrated Survival"
        )
      return(p)
    } else {
      p <- ggplot2::ggplot(pdp_integrated, ggplot2::aes_string(x = feature, y = "integrated_surv")) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste("Integrated PDP over Survival Time:", feature),
          x = feature,
          y = "Integrated Survival"
        )
      if (smooth) {
        p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "steelblue", size = 1.2)
      } else {
        p <- p + ggplot2::geom_line(color = "steelblue", size = 1.2)
      }
      return(p)
    }
  }

  # Per-time PDP/ICE plot (unchanged)
  plot_data <- results[type %in% switch(method,
                                        "pdp" = "pdp",
                                        "ice" = "ice",
                                        "pdp+ice" = c("pdp", "ice"),
                                        stop("Invalid method"))]

  if (!is.null(ids) && "ice" %in% plot_data$type) {
    plot_data <- plot_data[type != "ice" | .id %in% ids]
  }

  if (is_categorical) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = feature, y = "surv_prob", fill = "type")) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        fill = "Type"
      )

    if ("ice" %in% plot_data$type) {
      p <- p + ggplot2::geom_boxplot(data = plot_data[type == "ice"], alpha = alpha_ice, position = "dodge")
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + ggplot2::stat_summary(data = plot_data[type == "pdp"],
                                     fun = mean, geom = "point",
                                     shape = 21, fill = "black", size = 3,
                                     position = ggplot2::position_dodge(width = 0.75))
    }

  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = feature, y = "surv_prob",
                                                        group = "interaction(.id, time)", color = "type")) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        color = NULL
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))

    if ("ice" %in% plot_data$type) {
      p <- p + ggplot2::geom_line(data = plot_data[type == "ice"], alpha = alpha_ice)
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + ggplot2::geom_line(data = plot_data[type == "pdp"], size = 1.2)
    }
  }

  return(p)
}

library(survival)
data(veteran)
veteran$celltype <- as.factor(veteran$celltype)

mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)

# Numerical variable
pdp_ice_result1 <- compute_pdp(
  model = mod_ranger,
  data = veteran,
  feature = "age",
  times = c(100, 200, 300),
  method = "pdp+ice"
)

plot_pdp(pdp_ice_result1, feature = "age", which = "per_time")
plot_pdp(pdp_ice_result1, feature = "age", which = "integrated", smooth = TRUE)

# Categorical variable
pdp_ice_result2 <- compute_pdp(
  model = mod_ranger,
  data = veteran,
  feature = "celltype",
  times = c(100, 300),
  method = "pdp+ice"
)

plot_pdp(pdp_ice_result2, feature = "celltype", which = "per_time")



