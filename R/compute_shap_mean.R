# compute mean SHAP across a set of observations
compute_shap_mean <- function(model, newdata, baseline_data, times,
                              sample.size = 100, method = "meanabs") {
  stopifnot(nrow(newdata) >= 1)
  
  # compute SHAP for each row in newdata
  shap_list <- lapply(seq_len(nrow(newdata)), function(i) {
    compute_shap(
      model         = model,
      newdata       = newdata[i, , drop = FALSE],
      baseline_data = baseline_data,
      times         = times,
      sample.size   = sample.size,
      aggregate     = TRUE,
      method        = method
    )
  })
  
  # combine and average phi across rows
  shap_df <- data.table::rbindlist(shap_list, idcol = "id", fill = TRUE)
  shap_mean <- shap_df[, list(phi = mean(phi, na.rm = TRUE)), by = feature]
  
  attr(shap_mean, "shap_method") <- method
  return(shap_mean)
}
