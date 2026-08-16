.pmap_rbind_dt <- function(grid, f) {
  grid <- as.data.frame(grid)
  results <- lapply(seq_len(nrow(grid)), function(i) {
    do.call(f, as.list(grid[i, , drop = FALSE]))
  })
  data.table::rbindlist(results, fill = TRUE)
}

.complete_cases_df <- function(data, vars) {
  DT <- data.table::as.data.table(data)
  as.data.frame(DT[stats::complete.cases(DT[, vars, with = FALSE])])
}

.wide_metric_row <- function(params, cv_summary_dt) {
  metric_vals <- as.list(stats::setNames(cv_summary_dt$mean, cv_summary_dt$metric))
  do.call(data.table::data.table, c(params, metric_vals))
}

.arrange_by_metric_dt <- function(dt, metric, maximize) {
  dt <- data.table::as.data.table(dt)
  data.table::setorderv(dt, metric, order = if (maximize) -1L else 1L)
  dt[]
}
