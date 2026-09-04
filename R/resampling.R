#' Resampling specifications for evaluate(), tune() and compare()
#'
#' These constructors return a lightweight `survalis_resampling` description that
#' the verb-level functions turn into concrete train/test splits. They do not
#' touch the data themselves, so the same specification can be reused across
#' models to obtain paired comparisons (identical folds via a shared `seed`).
#'
#' * `cv()` -- k-fold cross-validation, optionally repeated, stratified on the
#'   event indicator by default.
#' * `holdout()` -- a single stratified train/test split.
#' * `group_cv()` -- k-fold cross-validation where whole groups (e.g. centres,
#'   subjects with repeated rows) are held out together.
#' * `bootstrap()` -- resampling with replacement; the out-of-bag rows form the
#'   assessment set (Efron-Tibshirani style, suitable for `.632` reporting).
#'
#' @param v Integer number of folds (default `5`).
#' @param repeats Integer number of times the whole v-fold scheme is repeated
#'   with fresh random folds (default `1`).
#' @param strata Column name to stratify on, or `NA` to disable stratification.
#'   The default sentinel `"status"` means "use the event-indicator column of the
#'   model formula", resolved by the caller.
#' @param prop Proportion of rows assigned to the training set in `holdout()`
#'   (default `0.75`).
#' @param group Column name whose distinct values define the groups held out
#'   together in `group_cv()`.
#' @param times Integer number of bootstrap resamples (default `25`).
#' @param seed Optional integer seed. Set it (and keep it equal across models)
#'   for paired resampling.
#'
#' @return An object of class `survalis_resampling`.
#'
#' @examples
#' cv(5)
#' cv(v = 10, repeats = 3, seed = 1)
#' holdout(prop = 0.8)
#' bootstrap(times = 50)
#' @name resampling
NULL

.new_resampling <- function(type, params) {
  structure(list(type = type, params = params), class = "survalis_resampling")
}

#' @rdname resampling
#' @export
cv <- function(v = 5, repeats = 1, strata = "status", seed = NULL) {
  v <- as.integer(v)
  repeats <- as.integer(repeats)
  if (v < 2L) stop("`v` must be >= 2.", call. = FALSE)
  if (repeats < 1L) stop("`repeats` must be >= 1.", call. = FALSE)
  .new_resampling("cv", list(v = v, repeats = repeats, strata = strata, seed = seed))
}

#' @rdname resampling
#' @export
holdout <- function(prop = 0.75, strata = "status", seed = NULL) {
  if (!is.finite(prop) || prop <= 0 || prop >= 1) {
    stop("`prop` must be in (0, 1).", call. = FALSE)
  }
  .new_resampling("holdout", list(prop = prop, strata = strata, seed = seed))
}

#' @rdname resampling
#' @export
group_cv <- function(group, v = 5, seed = NULL) {
  if (missing(group) || !is.character(group) || length(group) != 1L) {
    stop("`group` must be a single column name.", call. = FALSE)
  }
  v <- as.integer(v)
  if (v < 2L) stop("`v` must be >= 2.", call. = FALSE)
  .new_resampling("group_cv", list(group = group, v = v, seed = seed))
}

#' @rdname resampling
#' @export
bootstrap <- function(times = 25, seed = NULL) {
  times <- as.integer(times)
  if (times < 1L) stop("`times` must be >= 1.", call. = FALSE)
  .new_resampling("bootstrap", list(times = times, seed = seed))
}

#' @export
print.survalis_resampling <- function(x, ...) {
  p <- x$params
  desc <- switch(x$type,
    cv = sprintf("%d-fold CV%s", p$v,
                 if (p$repeats > 1L) sprintf(" x %d repeats", p$repeats) else ""),
    holdout = sprintf("holdout (%.0f%% train)", 100 * p$prop),
    group_cv = sprintf("%d-fold grouped CV on '%s'", p$v, p$group),
    bootstrap = sprintf("bootstrap (%d resamples, OOB assessment)", p$times),
    x$type
  )
  strata <- p$strata
  extra <- if (!is.null(strata) && !is.na(strata)) sprintf(", strata = %s", strata) else ""
  seed <- if (!is.null(p$seed)) sprintf(", seed = %s", p$seed) else ""
  cat(sprintf("<survalis_resampling> %s%s%s\n", desc, extra, seed))
  invisible(x)
}

#' Expand a resampling specification into concrete row-index splits
#'
#' @param resampling A `survalis_resampling` object.
#' @param data The (already row-filtered) modelling data frame.
#' @param status_col Name of the event-indicator column, used when
#'   `params$strata == "status"`.
#' @return A list of `list(analysis = <int>, assessment = <int>, id = <chr>)`.
#' @keywords internal
.make_splits <- function(resampling, data, status_col = NULL) {
  stopifnot(inherits(resampling, "survalis_resampling"))
  p <- resampling$params
  if (!is.null(p$seed)) set.seed(p$seed)

  strata_col <- NULL
  if (!is.null(p$strata) && !is.na(p$strata)) {
    strata_col <- if (identical(p$strata, "status")) status_col else p$strata
    if (!is.null(strata_col) && !strata_col %in% names(data)) strata_col <- NULL
  }

  to_idx <- function(split) {
    list(
      analysis = as.integer(split$in_id),
      assessment = as.integer(rsample::complement(split))
    )
  }

  rset <- switch(resampling$type,
    cv = if (is.null(strata_col)) {
      rsample::vfold_cv(data, v = p$v, repeats = p$repeats)
    } else {
      rsample::vfold_cv(data, v = p$v, repeats = p$repeats,
                        strata = !!rlang::sym(strata_col))
    },
    holdout = {
      sp <- if (is.null(strata_col)) {
        rsample::initial_split(data, prop = p$prop)
      } else {
        rsample::initial_split(data, prop = p$prop, strata = !!rlang::sym(strata_col))
      }
      rsample::manual_rset(list(sp), ids = "holdout")
    },
    group_cv = rsample::group_vfold_cv(data, group = !!rlang::sym(p$group), v = p$v),
    bootstrap = rsample::bootstraps(data, times = p$times),
    stop("Unknown resampling type: ", resampling$type)
  )

  ids <- if ("id2" %in% names(rset)) paste(rset$id, rset$id2, sep = ".") else rset$id
  Map(function(sp, id) c(to_idx(sp), list(id = id)), rset$splits, ids)
}
