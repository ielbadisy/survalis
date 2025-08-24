#' Veteran's Administration Lung Cancer Trial Data
#'
#' This is the `veteran` dataset originally from the \pkg{survival} package,
#' containing data from a randomized trial of lung cancer treatments.
#'
#' @format A data frame with 137 observations and 8 variables:
#' \describe{
#'   \item{trt}{Treatment: 1=standard, 2=test}
#'   \item{celltype}{Cell type: squamous, smallcell, adeno, large}
#'   \item{time}{Survival time in days}
#'   \item{status}{Censoring status: 1=dead, 0=alive}
#'   \item{karno}{Karnofsky performance score (higher = better)}
#'   \item{diagtime}{Months from diagnosis to randomization}
#'   \item{age}{Age in years}
#'   \item{prior}{Prior therapy: 0=no, 10=yes}
#' }
#'
#' @source \pkg{survival} package, originally from
#' Kalbfleisch and Prentice (1980) *The Statistical Analysis of Failure Time Data*.
"veteran"
