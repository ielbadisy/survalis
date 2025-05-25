fit_rpart <- function(formula, data,
                      minsplit = 20,
                      minbucket = round(minsplit / 3),
                      cp = 0.01,
                      maxcompete = 4,
                      maxsurrogate = 5,
                      usesurrogate = 2,
                      xval = 10,
                      surrogatestyle = 0,
                      maxdepth = 30) {
  stopifnot(requireNamespace("rpart", quietly = TRUE))

  method <- "exp"  # fixed to survival tree with exponential splitting

  model <- rpart::rpart(
    formula,
    data = data,
    method = method,
    control = list(
      minsplit = minsplit, minbucket = minbucket, cp = cp,
      maxcompete = maxcompete, maxsurrogate = maxsurrogate,
      usesurrogate = usesurrogate, xval = xval,
      surrogatestyle = surrogatestyle, maxdepth = maxdepth
    ),
    model = TRUE, x = TRUE, y = TRUE
  )

  structure(list(
    model = model,
    learner = "rpart",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2]
  ), class = "mlsurv_model", engine = "rpart")
}

predict_rpart <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "rpart")
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  model_party <- partykit::as.party(object$model)

  # predict expected survival times from the tree
  predicted_times <- stats::predict(model_party, newdata = newdata, type = "response")

  # convert predicted survival times to survival probabilities using exponential assumption
  surv_mat <- sapply(times, function(t) {
    exp(-t / predicted_times)
  })

  if (is.vector(surv_mat)) surv_mat <- matrix(surv_mat, ncol = length(times))

  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(newdata)))
  colnames(surv_mat) <- paste0("t=", times)

  return(as.data.frame(surv_mat))
}
library(survival)
library(rpart)
library(partykit)

mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)

pred_rpart <- predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300, 400, 500))

