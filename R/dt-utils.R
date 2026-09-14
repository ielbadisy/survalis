.pmap_rbind_dt <- function(grid, f) {
  grid <- as.data.frame(grid)
  results <- lapply(seq_len(nrow(grid)), function(i) {
    do.call(f, as.list(grid[i, , drop = FALSE]))
  })
  .rbind_fill_dt(results)
}

.map_rbind_dt <- function(x, f) {
  results <- lapply(x, f)
  .rbind_fill_dt(results)
}

# Row-binds `results` (data.frames, possibly with differing columns --
# e.g. a failed tuning row that omits metric columns -- and possibly
# carrying a list-column, e.g. survdnn's `hidden`) via base rbind() after
# padding to a common column set. Not basetable::rbindfill(): rbindfill()
# coerces list-columns to NA instead of preserving them.
.rbind_fill_dt <- function(results) {
  results <- Filter(Negate(is.null), results)
  if (!length(results)) return(data.frame())
  results <- lapply(results, as.data.frame)
  all_cols <- unique(unlist(lapply(results, names)))
  results <- lapply(results, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]
  })
  do.call(rbind, results)
}

.complete_cases_df <- function(data, vars) {
  data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
}

.wide_metric_row <- function(params, cv_summary_dt) {
  metric_vals <- as.list(stats::setNames(cv_summary_dt$mean, cv_summary_dt$metric))
  # list-valued params (e.g. survdnn's `hidden`) must stay a single list-column
  # entry rather than being spread across columns by base data.frame()'s
  # default list-argument handling.
  params <- lapply(params, function(x) if (is.list(x)) I(x) else x)
  do.call(data.frame, c(params, metric_vals, stringsAsFactors = FALSE))
}

.arrange_by_metric_dt <- function(dt, metric, maximize) {
  basetable::orderrows(dt, by = metric, decreasing = maximize)
}

# A plain data.frame/tibble/basetable's `[i, j, drop = FALSE]` always
# selects columns when `j` is a variable holding column names.
.select_cols <- function(x, cols) {
  x[, cols, drop = FALSE]
}
