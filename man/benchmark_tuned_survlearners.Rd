\name{benchmark_tuned_survlearners}
\alias{benchmark_tuned_survlearners}
\title{Benchmark Tuned Survival Learners with Nested Cross-Validation}
\usage{
benchmark_tuned_survlearners(
  formula,
  data,
  learners,
  times,
  metrics = c("cindex", "ibs"),
  outer_folds = 5,
  inner_folds = 5,
  seed = 123,
  inner_ncores = 1,
  learner_args = list(),
  refit_final = FALSE,
  verbose = FALSE,
  suppress_errors = TRUE,
  ...
)
}
\arguments{
\item{formula}{A survival formula of the form \code{Surv(time, status) ~ x1 + x2 + ...}.}

\item{data}{A data frame containing the variables in \code{formula}.}

\item{learners}{Character vector of learner ids (without prefixes), e.g. \code{c("ranger", "glmnet")}.}

\item{times}{Numeric vector of evaluation time points for survival predictions.}

\item{metrics}{Character vector of metrics to optimize and report.}

\item{outer_folds}{Integer number of outer CV folds used for performance estimation.}

\item{inner_folds}{Integer number of inner CV folds used by each learner's tuning routine.}

\item{seed}{Integer random seed used for outer and inner resampling.}

\item{inner_ncores}{Integer; number of CPU cores passed to each learner's inner
\code{cv_survlearner()} calls during tuning (default \code{1}).}

\item{learner_args}{Optional named list of learner-specific arguments. Each entry can be either a list of tuning arguments passed to \code{tune_<learner>()}, or a list with \code{tune} and/or \code{fit} components. The \code{fit} component is forwarded only when refitting the selected hyperparameters on an outer training split or on the full dataset.}

\item{refit_final}{Logical; if \code{TRUE}, tunes each learner on the full dataset and refits a final model with the selected hyperparameters.}

\item{verbose}{Logical; if \code{TRUE}, prints progress per learner and outer fold.}

\item{suppress_errors}{Logical; if \code{TRUE} (default), errors from individual learners or outer folds are caught and reported via \code{warning()}.}

\item{...}{Additional arguments passed to each learner's \code{tune_*()} function.}
}
\value{
A list of class \code{"nested_surv_benchmark"} with components \code{outer_results},
\code{outer_summary}, \code{selected_params}, \code{final_models}, and \code{settings}.
}
\description{
Runs nested cross-validation for one or more learners that expose
\code{tune_<learner>()}, \code{fit_<learner>()}, and \code{predict_<learner>()}.
The outer folds estimate performance, while each inner training split is tuned
using the learner's existing \code{tune_*()} implementation only on the outer
training data.
}
\examples{
\dontrun{
res <- benchmark_tuned_survlearners(
  Surv(time, status) ~ age + karno + trt,
  data = veteran,
  learners = c("ranger", "glmnet"),
  times = c(100, 200),
  outer_folds = 3,
  inner_folds = 2
)
res$outer_summary
}
}
\seealso{
\code{\link{benchmark_default_survlearners}}, \code{\link{cv_survlearner}}
}
