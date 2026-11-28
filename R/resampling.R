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
#' * `groupcv()` -- k-fold cross-validation where whole groups (e.g. centres,
#'   subjects with repeated rows) are held out together.
#' * `bootstrap()` -- resampling with replacement; the out-of-bag rows form the
#'   assessment set. [evaluate()] reports the out-of-bag estimate or, with
#'   `estimator`, the .632, .632+ or optimism-corrected estimate (see Details).
#'
#' @details
#' Bootstrap estimators. Let \eqn{\overline{err}} be the apparent score (the
#' learner fitted and scored on all the data), \eqn{Err^{(1)}} the out-of-bag
#' score (mean over the resamples of the score of fit \eqn{b} on the rows left out
#' of resample \eqn{b}) and \eqn{\gamma} the no-information score, the score of the
#' full-data predictions paired with outcomes at random, averaged over all
#' \eqn{n^2} pairs (Efron and Tibshirani 1997, eq. 27). For the Brier score and the
#' integrated Brier score \eqn{\gamma} is computed exactly at each time (the score
#' is quadratic in the prediction, with IPCW weights that follow the outcome, as in
#' Gerds and Schumacher 2007); for the other metrics it is the mean score over
#' `nperm` random permutations of the prediction rows, an unbiased estimate of
#' the \eqn{n^2}-pair mean. Then
#' * `".632"`: \eqn{0.368\,\overline{err} + 0.632\,Err^{(1)}}.
#' * `".632+"`: \eqn{(1 - w)\,\overline{err} + w\,Err^{(1)}} with
#'   \eqn{w = 0.632 / (1 - 0.368 R)} and relative overfitting rate
#'   \eqn{R = (Err^{(1)} - \overline{err}) / (\gamma - \overline{err})}, set to 0 when
#'   \eqn{Err^{(1)} \le \overline{err}} or \eqn{\gamma \le \overline{err}} and to 1
#'   when \eqn{Err^{(1)} \ge \gamma} (Efron and Tibshirani 1997, eqs. 28-32; the
#'   censored-data version of Gerds and Schumacher 2007, Section 4). With
#'   \eqn{R = 1} the estimate is \eqn{Err^{(1)}}.
#' * `"optimism"`: \eqn{\overline{err}} plus the mean over resamples of the score
#'   of fit \eqn{b} on all the data minus its score on its own resample (Harrell's
#'   optimism bootstrap).
#'
#' Each estimator is applied to each metric separately (per time for
#' `"brier@<time>"`, to the integrated score for `"ibs"`). For a metric where
#' higher is better (`"cindex"`, `"auc"`) the formulas are applied to its negative,
#' so the inequalities are reversed. The same bootstrap resamples serve all four
#' estimators, and [evaluate()] returns all of them in its `bootstrap` element.
#' [tune()] and [compare()] score the out-of-bag rows only.
#'
#' References: Efron B, Tibshirani R (1997). Improvements on cross-validation: the
#' .632+ bootstrap method. JASA 92, 548-560. Gerds TA, Schumacher M (2007).
#' Efron-type measures of prediction error for survival analysis. Biometrics 63,
#' 1283-1287.
#'
#' @param v Integer number of folds (default `5`).
#' @param repeats Integer number of times the whole v-fold scheme is repeated
#'   with fresh random folds (default `1`).
#' @param strata Column name to stratify on, or `NA` to disable stratification.
#'   The default sentinel `"status"` means "use the event-indicator column of the
#'   model formula", whatever its name (for example `cens` in
#'   `Surv(time, cens) ~ .`). Any other name must be a column of the data;
#'   otherwise a warning is given and the resampling is not stratified.
#' @param prop Proportion of rows assigned to the training set in `holdout()`
#'   (default `0.75`).
#' @param group Column name whose distinct values define the groups held out
#'   together in `groupcv()`.
#' @param times Integer number of bootstrap resamples (default `25`).
#' @param estimator Bootstrap estimate that [evaluate()] reports as the `mean` of
#'   its summary: `"oob"` (out-of-bag, default), `".632"`, `".632+"` or
#'   `"optimism"`. Any choice other than `"oob"` costs one extra fit on all the data
#'   and a prediction of all the rows by each bootstrap fit.
#' @param nperm Number of random permutations for the no-information score of the
#'   metrics other than the Brier scores (default `20`).
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
#' bootstrap(times = 50, estimator = ".632+")
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
groupcv <- function(group, v = 5, seed = NULL) {
  if (missing(group) || !is.character(group) || length(group) != 1L) {
    stop("`group` must be a single column name.", call. = FALSE)
  }
  v <- as.integer(v)
  if (v < 2L) stop("`v` must be >= 2.", call. = FALSE)
  .new_resampling("groupcv", list(group = group, v = v, seed = seed))
}

#' @rdname resampling
#' @export
bootstrap <- function(times = 25, seed = NULL, estimator = c("oob", ".632", ".632+", "optimism"),
                      nperm = 20) {
  times <- as.integer(times)
  if (times < 1L) stop("`times` must be >= 1.", call. = FALSE)
  .new_resampling("bootstrap", list(times = times, seed = seed, estimator = match.arg(estimator),
                                    nperm = as.integer(nperm)))
}

#' @export
print.survalis_resampling <- function(x, ...) {
  p <- x$params
  desc <- switch(x$type,
    cv = sprintf("%d-fold CV%s", p$v,
                 if (p$repeats > 1L) sprintf(" x %d repeats", p$repeats) else ""),
    holdout = sprintf("holdout (%.0f%% train)", 100 * p$prop),
    groupcv = sprintf("%d-fold grouped CV on '%s'", p$v, p$group),
    bootstrap = sprintf("bootstrap (%d resamples, %s estimate)", p$times,
                        .null_default(p$estimator, "oob")),
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
    if (!is.null(strata_col) && !strata_col %in% names(data)) {
      warning("Stratification column ", sQuote(strata_col), " is not in the data; ",
              "resampling is not stratified. Use `strata = \"status\"` (the default) ",
              "to stratify on the event column of the formula, or `strata = NA`.",
              call. = FALSE)
      strata_col <- NULL
    }
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
    groupcv = rsample::group_vfold_cv(data, group = !!rlang::sym(p$group), v = p$v),
    bootstrap = rsample::bootstraps(data, times = p$times),
    stop("Unknown resampling type: ", resampling$type)
  )

  ids <- if ("id2" %in% names(rset)) paste(rset$id, rset$id2, sep = ".") else rset$id
  Map(function(sp, id) c(to_idx(sp), list(id = id)), rset$splits, ids)
}
