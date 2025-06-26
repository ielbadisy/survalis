# New base R + dplyr version of compute_shap
compute_shap <- function(model, predict_function, newdata, baseline_data,
                         times, sample.size = 100,
                         aggregate = FALSE, method = c("meanabs", "integral")) {
  stopifnot(nrow(newdata) == 1)
  method <- match.arg(method)

  feature_names <- setdiff(colnames(newdata), c("time", "status", "event"))
  n_features <- length(feature_names)
  all_results <- list()

  for (current_time in times) {
    shapley_values <- replicate(sample.size, {
      order_features <- sample(feature_names)
      baseline_sample <- baseline_data[sample(nrow(baseline_data), 1), , drop = FALSE]
      instance_with <- newdata
      instance_without <- newdata
      contributions <- numeric(n_features)

      for (k in seq_len(n_features)) {
        feature_k <- order_features[k]
        if (k < n_features) {
          future_features <- order_features[(k + 1):n_features]
          instance_with[future_features] <- baseline_sample[future_features]
        }
        instance_without <- instance_with
        instance_without[feature_k] <- baseline_sample[[feature_k]]

        pred_with <- predict_function(model, instance_with, times)
        pred_without <- predict_function(model, instance_without, times)

        idx_time <- which.min(abs(times - current_time))
        pred_with <- if (is.vector(pred_with)) pred_with[idx_time] else pred_with[, idx_time]
        pred_without <- if (is.vector(pred_without)) pred_without[idx_time] else pred_without[, idx_time]

        idx_feature <- which(feature_names == feature_k)
        contributions[idx_feature] <- contributions[idx_feature] + (pred_with - pred_without)
      }
      contributions
    }, simplify = TRUE)

    shapley_df <- as.data.frame(t(shapley_values))
    colnames(shapley_df) <- feature_names
    mean_phi <- colMeans(shapley_df)
    var_phi <- apply(shapley_df, 2, var)

    result_df <- data.frame(
      feature = feature_names,
      phi = as.numeric(mean_phi),
      phi.var = as.numeric(var_phi),
      time = current_time
    )
    all_results[[as.character(current_time)]] <- result_df
  }

  final_result <- do.call(rbind, all_results)

  if (aggregate) {
    library(tidyr)
    library(dplyr)

    phi_wide <- final_result %>%
      select(feature, time, phi) %>%
      pivot_wider(names_from = time, values_from = phi)

    phi_matrix <- as.matrix(phi_wide[, -1])
    rownames(phi_matrix) <- phi_wide$feature

    if (method == "meanabs") {
      agg_phi <- apply(phi_matrix, 1, function(row) mean(abs(row)))
    } else {
      agg_phi <- apply(phi_matrix, 1, function(row) pracma::trapz(times, row))
    }

    phi_var <- apply(phi_matrix, 1, var)
    result_agg <- data.frame(
      feature = rownames(phi_matrix),
      phi = agg_phi,
      phi.var = phi_var
    )
    attr(result_agg, "shap_method") <- method
    return(result_agg)
  }

  return(final_result)
}



plot_shap <- function(shapley_result, type = c("auto", "time", "integrated")) {
  type <- match.arg(type)

  is_time_dependent <- "time" %in% colnames(shapley_result)

  if (type == "auto") {
    type <- if (is_time_dependent) "time" else "integrated"
  }

  if (type == "time") {
    if (!is_time_dependent) stop("No time column found in SHAP object. Set type = 'integrated'.")

    ggplot2::ggplot(shapley_result, ggplot2::aes(x = time, y = phi, color = feature)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_ribbon(
        data = shapley_result[!is.na(shapley_result$phi.var), ],
        ggplot2::aes(ymin = phi - sqrt(phi.var), ymax = phi + sqrt(phi.var), fill = feature),
        alpha = 0.2,
        color = NA,
        inherit.aes = TRUE
      ) +
      ggplot2::labs(
        title = "Shapley Contributions Over Time",
        x = "Time",
        y = "SHAP value"
      ) +
      ggplot2::theme_minimal(base_size = 14)
  }

  if (type == "integrated") {
    shap_method <- attr(shapley_result, "shap_method", exact = TRUE)
    shap_method <- if (is.null(shap_method)) "integral" else shap_method

    label <- switch(
      shap_method,
      "meanabs" = "Mean Absolute SHAP",
      "integral" = "Integrated SHAP",
      "SHAP Summary"
    )

    shapley_result <- shapley_result[order(shapley_result$phi), ]
    shapley_result$feature <- factor(shapley_result$feature, levels = shapley_result$feature)

    p <- ggplot2::ggplot(shapley_result, ggplot2::aes(x = phi, y = feature)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(
        title = label,
        x = "SHAP value",
        y = "Feature"
      ) +
      ggplot2::theme_minimal(base_size = 14)

    # Add error bars if variance exists
    if ("phi.var" %in% names(shapley_result)) {
      p <- p + ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = phi - sqrt(phi.var), xmax = phi + sqrt(phi.var)),
        height = 0.2
      )
    }

    return(p)
  }
}


library(survival)
library(ranger)
library(pracma)
library(reshape2)

data(veteran)
times <- c(100, 200, 300)
obs <- veteran[100, , drop = FALSE]

mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran, num.trees = 300)

shap_agg_meanabs <- compute_shap(
  model = mod_ranger,
  predict_function = predict_ranger,
  newdata = obs,
  baseline_data = veteran,
  times = times,
  sample.size = 30,
  aggregate = TRUE,
  method = "meanabs"
)

shap_agg_integral <- compute_shap(
  model = mod_ranger,
  predict_function = predict_ranger,
  newdata = obs,
  baseline_data = veteran,
  times = times,
  sample.size = 30,
  aggregate = TRUE,
  method = "integral"
)




plot_shap(shap_agg_meanabs)
plot_shap(shap_agg_integral)

