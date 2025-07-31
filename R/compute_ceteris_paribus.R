compute_ceteris_paribus <- function(model, predict_function, data, instances, feature, times, grid.size = 20) {
  stopifnot(feature %in% names(data))

  is_categorical <- is.factor(data[[feature]]) || is.character(data[[feature]])
  grid_values <- if (is_categorical) {
    sort(unique(data[[feature]]))
  } else {
    seq(min(data[[feature]], na.rm = TRUE), max(data[[feature]], na.rm = TRUE), length.out = grid.size)
  }

  all_profiles <- purrr::map_df(seq_len(nrow(instances)), function(i) {
    instance <- instances[i, , drop = FALSE]
    replicated_data <- do.call(rbind, replicate(length(grid_values), instance, simplify = FALSE))
    replicated_data[[feature]] <- grid_values

    preds <- predict_function(model, newdata = replicated_data, times = times)
    if (is.vector(preds)) preds <- matrix(preds, ncol = 1)

    df <- as.data.frame(preds)
    colnames(df) <- paste0("t=", times)
    df[[feature]] <- grid_values
    df$.id <- i

    reshape(df,
            varying = paste0("t=", times),
            v.names = "surv_prob",
            timevar = "time",
            times = times,
            idvar = c(".id", feature),
            direction = "long")
  })

  all_profiles$type <- "ceteris_paribus"
  rownames(all_profiles) <- NULL
  return(all_profiles)
}

plot_ceteris_paribus <- function(cp_data, feature, alpha = 0.4, facet = FALSE) {
  library(ggplot2)
  is_categorical <- is.factor(cp_data[[feature]]) || is.character(cp_data[[feature]])

  p <- ggplot(cp_data, aes_string(x = feature, y = "surv_prob", group = "interaction(.id, time)", color = "as.factor(time)")) +
    geom_line(alpha = alpha) +
    labs(
      x = feature,
      y = "Survival Probability",
      title = paste("Ceteris Paribus Profiles for", feature),
      color = "Time"
    ) +
    theme_minimal()

  if (facet) {
    p <- p + facet_wrap(~ time)
  }

  return(p)
}

# select multiple individuals
subset_rows <- veteran[1:200, ]

cp_result <- compute_ceteris_paribus(
  model = mod_cox,
  predict_function = predict_cox,
  data = veteran,
  instances = subset_rows,
  feature = "age",
  times = c(100, 200, 300, 400, 500)
)

plot_ceteris_paribus(cp_result, feature = "age")

plot_ceteris_paribus(cp_result, feature = "age", facet = TRUE)
