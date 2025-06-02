compute_shap <- function(model, predict_function, newdata, baseline_data,
                         times, sample.size = 100,
                         aggregate = FALSE, method = c("meanabs", "integral")) {
  requireNamespace("data.table")
  requireNamespace("pracma")
  
  method <- match.arg(method)
  
  if (nrow(newdata) != 1) stop("newdata must contain exactly one observation.")
  
  feature_names <- setdiff(colnames(newdata), c(model$time, model$status))
  n_features <- length(feature_names)
  
  all_results <- list()
  
  for (current_time in times) {
    ## use multicore here to reduce computing time 
    shapley_values <- lapply(seq_len(sample.size), function(i) {
      order_features <- sample(feature_names)
      baseline_sample <- baseline_data[sample(nrow(baseline_data), 1), , drop = FALSE]
      
      instance_with <- newdata
      instance_without <- newdata
      contributions <- numeric(n_features)
      
      for (k in seq_len(n_features)) {
        feature_k <- order_features[k]
        
        if (k < n_features) {
          future_features <- order_features[(k + 1):n_features]
          instance_with[, future_features] <- baseline_sample[, future_features]
        }
        
        instance_without <- instance_with
        instance_without[, feature_k] <- baseline_sample[, feature_k]
        
        pred_with <- predict_function(model, instance_with, times = times)
        pred_without <- predict_function(model, instance_without, times = times)
        
        if (is.null(dim(pred_with))) pred_with <- matrix(pred_with, nrow = 1)
        if (is.null(dim(pred_without))) pred_without <- matrix(pred_without, nrow = 1)
        
        idx_time <- which.min(abs(times - current_time))
        
        pred_with <- pred_with[, idx_time]
        pred_without <- pred_without[, idx_time]
        
        idx_feature <- which(feature_names == feature_k)
        contributions[idx_feature] <- contributions[idx_feature] + (pred_with - pred_without)
      }
      
      contributions
    })
    
    shapley_dt <- data.table::as.data.table(do.call(rbind, shapley_values))
    data.table::setnames(shapley_dt, feature_names)
    
    mean_dt <- shapley_dt[, lapply(.SD, mean)]
    var_dt <- shapley_dt[, lapply(.SD, var)]
    
    result_dt <- data.table::data.table(
      feature = feature_names,
      phi = as.numeric(mean_dt),
      phi.var = as.numeric(var_dt),
      time = current_time
    )
    
    all_results[[as.character(current_time)]] <- result_dt
  }
  
  final_result <- data.table::rbindlist(all_results)
  
  # optional aggregation across time
  if (aggregate) {
    phi_wide <- reshape2::dcast(final_result, feature ~ time, value.var = "phi")
    rownames(phi_wide) <- phi_wide$feature
    phi_wide$feature <- NULL
    
    if (method == "meanabs") {
      agg_phi <- apply(phi_wide, 1, function(row) mean(abs(row)))
    } else if (method == "integral") {
      agg_phi <- apply(phi_wide, 1, function(row) pracma::trapz(times, row))
    }
    
    aggregated_result <- data.table::data.table(
      feature = rownames(phi_wide),
      phi = agg_phi
    )
    return(aggregated_result)
  }
  
  return(final_result)
}


library(survival)
library(ranger)
library(pracma)
library(reshape2)

data(veteran)
times <- c(100, 200, 300)
obs <- veteran[100, , drop = FALSE]

mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran, num.trees = 300)

# aggregated SHAP using mean absolute value
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

# aggregated SHAP using integral
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

print(shap_agg_meanabs)
print(shap_agg_integral)


plot_shap <- function(shapley_result, type = c("auto", "time", "integrated")) {
  type <- match.arg(type)
  
  # auto-detect if "time" column exists (time-dependent SHAP)
  is_time_dependent <- "time" %in% colnames(shapley_result)
  
  if (type == "auto") {
    type <- if (is_time_dependent) "time" else "integrated"
  }
  
  if (type == "time") {
    if (!is_time_dependent) stop("No time column found in SHAP object. Set type = 'integrated'.")
    
    ggplot2::ggplot(shapley_result, ggplot2::aes(x = time, y = phi, color = feature)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = phi - sqrt(phi.var), ymax = phi + sqrt(phi.var)),
                           alpha = 0.2, color = NA) +
      ggplot2::labs(
        title = "Shapley Contributions Over Time",
        x = "Time",
        y = "SHAP value"
      ) +
      ggplot2::theme_minimal(base_size = 14)
    
  } else if (type == "integrated") {
    if (!is_time_dependent) {
      # Single-time SHAP already aggregated (e.g., phi only)
      shapley_result <- shapley_result[order(shapley_result$phi), ]
      shapley_result$feature <- factor(shapley_result$feature, levels = shapley_result$feature)
      
      return(
        ggplot2::ggplot(shapley_result, ggplot2::aes(x = phi, y = feature)) +
          ggplot2::geom_col(fill = "steelblue") +
          ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = phi - sqrt(phi.var), xmax = phi + sqrt(phi.var)),
            height = 0.2
          ) +
          ggplot2::labs(
            title = "Shapley Contributions (Single Time)",
            x = "SHAP value",
            y = "Feature"
          ) +
          ggplot2::theme_minimal(base_size = 14)
      )
    }
    
    # Integrate SHAP(t) over time
    if (!requireNamespace("pracma", quietly = TRUE)) stop("Package 'pracma' required. Install via install.packages('pracma')")
    library(dplyr)
    
    integrated <- shapley_result |>
      group_by(feature) |>
      summarise(
        integrated_phi = pracma::trapz(time, phi)
      ) |>
      arrange(desc(integrated_phi))
    
    ggplot2::ggplot(integrated, ggplot2::aes(x = integrated_phi, y = reorder(feature, integrated_phi))) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(
        title = "Integrated SHAP Contributions",
        x = "Integrated SHAP value",
        y = "Feature"
      ) +
      ggplot2::theme_minimal(base_size = 14)
  }
}

plot_shap(shap)                     # Auto mode
plot_shap(shap, type = "time")      # Time-varying SHAP
plot_shap(shap, type = "integrated")# Area-under-curve SHAP
