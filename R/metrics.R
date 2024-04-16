#' Compute averaged error of predicted target variable.
#'
#' @param family family of regression: "gaussian"
#' @param actual matrix of actual target variable
#' @param pred matrix of predicted target variable
#'
#' @export cv.trlasso.mse
#' @export
cv.trlasso.mse <- function(family, actual, pred) {
  if (family == "gaussian") {
    mse <- apply(pred, 2, function(v) {
      mean((v - actual)^2)
    })
  }
  return(mse)
}

#' Compute errors of predicted target variable.
#'
#' @param family family of regression: "gaussian"
#' @param actual_all matrix of actual target variable
#' @param pred_all matrix of predicted target variable
#' @param direct_prediction either corrected cross validation is used or not
#'
#' @export cv.trlasso.resid_all
#' @export
cv.trlasso.resid_all <- function(family, actual_all, pred_all, direct_prediction = FALSE) {
  if (family == "gaussian") {
    resid_all <- apply(pred_all, 2, function(v) {
      (actual_all - v)^2
    })
  }
  return(resid_all)
}
