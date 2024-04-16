#' Predict responses
#'
#' @param object trlasso model
#' @param newx matrix of explanatory variables
#' @param s selected lambda (default: all)
#' @param type of prediction for logistic regression: "response"(default) or "class"
#' @param ... parameters of predict function
#'
#' @examples
#' X <- as.matrix(iris[, 1:3])
#' y <- iris[, 4]
#' cv_fit <- cv.trlasso(X, y, nlambda=50, lambda.min.ratio=1e-2)
#' predict(cv_fit$fit, X)
#'
#' @export predict.trlasso
#' @export
predict.trlasso <- function(object, newx, s=NULL, type="response", ...) {
  varnames <- colnames(newx)
  samplenames <- rownames(newx)
  
  if (object$family=="gaussian") {
    if (is.null(s)) {
      t(apply(newx %*% object$beta, 1, function(v){v + as.numeric(object$a0)}))
    } else {
      newx %*% object$beta[, which(object$lambda==s)[1], drop=FALSE] + object$a0[which(object$lambda==s)[1]]
    }
  } else {
    stop("specified family is not supported")
  }
}