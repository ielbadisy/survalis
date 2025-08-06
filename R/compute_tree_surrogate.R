compute_tree_surrogate <- function(model, data, times, minsplit = 10, cp = 0.01) {
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Please install 'rpart'.")
  if (!requireNamespace("partykit", quietly = TRUE)) stop("Please install 'partykit'.")

  if (length(times) < 1) stop("At least one evaluation time must be specified.")

  # infer prediction function from model
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  is_dynamic <- length(times) > 1

  results <- lapply(times, function(time) {
    preds <- predict_function(model, newdata = data, times = time)

    # Handle multiple formats (data.frame, matrix, etc.)
    col_match <- grep(paste0("^t=*", time, "$"), colnames(preds))
    if (length(col_match) != 1) {
      warning(sprintf("Could not find a single matching column for time = %s. Using first column.", time))
      y <- if (is.data.frame(preds)) preds[[1]] else preds[, 1]
    } else {
      y <- preds[[col_match]]
    }

    # Exclude common survival columns
    exclude_vars <- c("time", "status", "event")
    covariates <- setdiff(names(data), exclude_vars)
    dat <- cbind(data[, covariates, drop = FALSE], surv_prob = y)

    # Fit surrogate tree
    tree_fit <- rpart::rpart(
      surv_prob ~ .,
      data = dat,
      control = rpart::rpart.control(minsplit = minsplit, cp = cp)
    )

    yhat <- predict(tree_fit, newdata = dat)
    sst <- var(y)
    sse <- var(y - yhat)
    r2 <- 1 - sse / sst

    split_count <- table(tree_fit$frame$var[tree_fit$frame$var != "<leaf>"])

    list(
      time = time,
      tree = tree_fit,
      r_squared = r2,
      split_count = split_count
    )
  })

  names(results) <- paste0("time_", times)

  structure(list(
    times = times,
    results = results,
    dynamic = is_dynamic
  ), class = "tree_surrogate")
}


mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)

tree_ranger <- compute_tree_surrogate(
  model = mod_ranger,
  data = veteran,
  times = c(100, 200, 300)
)

plot_tree_surrogate(tree_ranger, type = "tree")
plot_tree_surrogate(tree_ranger, type = "importance")


