compute_cate_rmst <- function(model, data, times, predict_fun) {
  data1 <- data; data1$A <- 1
  data0 <- data; data0$A <- 0

  S1 <- predict_fun(model, newdata = data1, times = times)
  S0 <- predict_fun(model, newdata = data0, times = times)

  rmst1 <- compute_rmst(S1, times)
  rmst0 <- compute_rmst(S0, times)

  data.frame(
    cate = rmst1 - rmst0,
    rmst_treated = rmst1,
    rmst_control = rmst0
  )
}


plot_cate <- function(cate_df, stratify_by = NULL, bins = 40) {
  if (!"cate" %in% names(cate_df)) {
    stop("`cate_df` must include a 'cate' column.")
  }

  library(ggplot2)

  if (!is.null(stratify_by)) {
    if (!stratify_by %in% names(cate_df)) {
      stop("stratify_by must be a column in cate_df.")
    }

    ggplot(cate_df, aes(x = cate, fill = factor(.data[[stratify_by]]))) +
      geom_density(alpha = 0.5) +
      labs(
        title = "CATE Distribution by Subgroup",
        x = "CATE (RMST difference)", y = "Density", fill = stratify_by
      ) +
      theme_minimal()
  } else {
    ggplot(cate_df, aes(x = cate)) +
      geom_histogram(bins = bins, fill = "darkblue", alpha = 0.7) +
      labs(
        title = "Distribution of Individual Treatment Effects (CATE)",
        x = "CATE (RMST difference)", y = "Count"
      ) +
      theme_minimal()
  }
}


# Compute individual-level CATEs
cate_df <- compute_cate_rmst(model, data = df, times = times, predict_fun = predict_xgboost)

# Option 1: Plot full distribution
plot_cate(cate_df)

# Option 2: Stratified by X2
cate_df$X2 <- df$X2
plot_cate(cate_df, stratify_by = "X2")
