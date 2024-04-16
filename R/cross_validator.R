#' Fit a model using a design matrix with cross validation
#'
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param beta_tilde vector of initial parameters
#' @param family family of regression: "gaussian" is only supported
#' @param nfolds the number of folds (ignored if foldid is specified)
#' @param foldid vector indicating id of fold for each sample
#' @param lambda.min.ratio ratio of max lambda and min lambda (ignored if lambda is specified)
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param alpha alpha
#' @param unit unit for cross validation error: "sample" (default) or "fold"
#' @param seed random seed of cross validation
#' @param eps convergence threshold (default 1e-4)
#' @param verbose whether output verbose warnings and messages (default FALSE)
#' @param standardize standardize or not (experimental, not safe)
#' @param penalty.factor.lasso factors applied to each coefficient in lasso penalty
#' @param penalty.factor.trlasso factors applied to each coefficient in trlasso penalty
#' @param penalty.factor.thresh lower limit for penalty factors (default 1e-2)
#' @param ... parameters of trlasso function
#'
#' @return lasso model
#' \item{fit}{lasso model with hole data}
#' \item{lambda.min}{lambda with minimum cross validation error}
#' \item{lambda.min.index}{index of lambda.min}
#' \item{lambda.1se}{largest lambda such that error is within 1 standard error of the minimum}
#' \item{lambda.1se.index}{index of lambda.1se}
#' \item{foldid}{fold id}
#' \item{cve}{cross validation error}
#' \item{cvse}{cross validation standard error}
#' \item{cvup}{cross validation error + standard error}
#' \item{cvlo}{cross validation error - standard error}
#' \item{fit.preval}{array of prevalidated fits}
#'
#' @examples
#' X <- as.matrix(iris[, 1:3])
#' y <- iris[, 4]
#' cv_fit <- cv.trlasso(X, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(cv_fit)
#' plot(cv_fit$fit)
#'
#' @export
cv.trlasso <- function(X, y, beta_tilde = NULL, family="gaussian",
                       nfolds = 10,
                       lambda.min.ratio = 1e-2,
                       nlambda = 100,
                       lambda = NULL,
                       alpha = 1,
                       foldid = NULL,
                       unit = "sample",
                       seed = NULL,
                       eps = 1e-4,
                       verbose = FALSE,
                       standardize = TRUE,
                       penalty.factor.lasso = NULL,
                       penalty.factor.trlasso = NULL,
                       penalty.factor.thresh = 1e-2,
                       ...) {
  
  if (verbose) {
    message("preprocess")
  }
  
  if (is.null(foldid)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    foldid <- sample(rep(1:nfolds, floor(nrow(X) / nfolds)), replace=FALSE)
    # foldid <- split(1:nrow(X), folds)
  } else {
    nfolds <- length(unique(foldid))
  }
  
  # constract standardized X and y
  n <- nrow(X)
  X_mean <- apply(X, 2, function(v){mean(v, na.rm=TRUE)})
  if (standardize) {
    X_sd <- apply(X, 2, function(v){
      sqrt((sum(!is.na(v))-1) / sum(!is.na(v)) * var(v, na.rm=TRUE))
    })
    X_sd[X_sd < 1e-4] <- 1e-4
  } else {
    X_sd <- rep(1, length=ncol(X))
  }
  X_tilde <- sweep(sweep(X, 2, X_mean, "-"), 2, X_sd, "/")
  X_tilde[is.na(X_tilde) | is.nan(X_tilde)] <- 1
  y_tilde <- y - mean(y, na.rm=TRUE)
  
  # construct covariance matrices
  gamma <- apply(X_tilde, 2, function(v){
    fill <- which(!is.na(v) & !is.na(as.vector(y_tilde)))
    g <- sum(v[fill] * y_tilde[fill]) / length(fill)
    if (is.nan(g) | is.na(g)) {g <- 0}
    return(g)
  })
  # gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
  
  # initialize penalty
  if (is.null(penalty.factor.lasso)) {
    penalty.factor.lasso <- rep(1, length=ncol(X))
  } else {
    penalty.factor.lasso[penalty.factor.lasso < penalty.factor.thresh] <- penalty.factor.thresh
  }
  if (is.null(penalty.factor.trlasso)) {
    penalty.factor.trlasso <- rep(1, length=ncol(X))
  } else {
    penalty.factor.trlasso[penalty.factor.trlasso < penalty.factor.thresh] <- penalty.factor.thresh
  }

  # construct a lambda sequence
  if (is.null(lambda)) {
    if (is.null(beta_tilde)) {
      beta_tilde_standard <- rep(0, length=ncol(X))
    } else {
      beta_tilde_standard <- beta_tilde * X_sd
    }

    idx_gamma <- (alpha*(penalty.factor.lasso+penalty.factor.trlasso)-penalty.factor.trlasso)>0
    gamma_l1 <- ifelse(sum(beta_tilde_standard==0 & idx_gamma)!=0,
                       max(abs(gamma)[beta_tilde_standard==0 & idx_gamma]),
                       0)
    gamma_l2 <- ifelse(sum(beta_tilde_standard>0 & idx_gamma)!=0,
                        max((-gamma / (alpha*(penalty.factor.lasso-penalty.factor.trlasso)+penalty.factor.trlasso))[beta_tilde_standard > 0 & idx_gamma],
                            (gamma /(alpha*(penalty.factor.lasso+penalty.factor.trlasso)-penalty.factor.trlasso))[beta_tilde_standard > 0 & idx_gamma]),
                        0)
    gamma_l3 <- ifelse(sum(beta_tilde_standard<0 & idx_gamma)!=0,
                         max((-gamma / (penalty.factor.trlasso - alpha*(penalty.factor.lasso+penalty.factor.trlasso)))[beta_tilde_standard < 0 & idx_gamma],
                             (gamma / (penalty.factor.trlasso + alpha*(penalty.factor.lasso-penalty.factor.trlasso)))[beta_tilde_standard < 0 & idx_gamma]),
                         0)

    r <- t(X_tilde) %*% (y_tilde - X_tilde %*% as.matrix(beta_tilde_standard, ncol=1)) / n
    idx_r <- (alpha*(penalty.factor.lasso+penalty.factor.trlasso)-penalty.factor.trlasso)<0
    r_l1 <-ifelse(sum(beta_tilde_standard==0 & idx_r)!=0,
                       max(abs(r)[beta_tilde_standard==0 & idx_r]),
                       0)
    r_l2 <- ifelse(sum(beta_tilde_standard>0 & idx_r)!=0,
                     max((-r / (penalty.factor.trlasso- alpha*(penalty.factor.lasso+penalty.factor.trlasso)))[beta_tilde_standard > 0 & idx_r],
                         (r / (penalty.factor.trlasso + alpha*(penalty.factor.lasso-penalty.factor.trlasso)))[beta_tilde_standard > 0 & idx_r]),
                     0)
    r_l3 <- ifelse(sum(beta_tilde_standard<0 & idx_r)!=0,
                     max((-r / (penalty.factor.trlasso + alpha*(penalty.factor.lasso-penalty.factor.trlasso)))[beta_tilde_standard < 0 & idx_r],
                         (r / (penalty.factor.trlasso - alpha*(penalty.factor.lasso+penalty.factor.trlasso)))[beta_tilde_standard < 0 & idx_r]),
                     0)

    idx_l <- (alpha*(penalty.factor.lasso+penalty.factor.trlasso)-penalty.factor.trlasso == 0)
    l1 <- ifelse(sum(idx_l)!=0, max(abs(gamma[idx_l])), 0)
    l2 <- ifelse(sum(idx_l)!=0, max(abs(r[idx_l])), 0)
    l_1_2 <- min(l1, l2)

    lambda_max <- max(gamma_l1, gamma_l2, gamma_l3, r_l1, r_l2, r_l3, l_1_2)

    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  }
  
  # cross validation
  if (verbose) {
    message("cv")
  }
  cv_fit <- list()
  for (i in 1:nfolds){
  # cv_fit <- purrr::map(1:nfolds, function(i) {
    if (verbose) {
      message("library")
    }
    #library(trlasso)
    if (verbose) {
      message(paste0("fold #", i))
    }
    
    train_id <- which(foldid != i)
    test_id <- which(foldid == i)
    
    # fit the model
    if (verbose) {
      message("start fold trlasso")
    }
    fit <- trlasso(X[train_id, ], y[train_id], beta_tilde=beta_tilde, family=family, lambda=lambda, alpha=alpha, 
                   penalty.factor.lasso=penalty.factor.lasso, penalty.factor.trlasso=penalty.factor.trlasso, 
                   penalty.factor.thresh=penalty.factor.thresh, eps=eps, verbose=verbose, ...)
    
    # evaluate the model
    pred <- predict.trlasso(fit, X[test_id, ])
    actual <- y[test_id]
    resid <- actual - pred
    # BUG: when actual & pred is vector (ex: test case is N=1)
    mse <- cv.trlasso.mse(family, actual, pred)
    
    cv_fit[[i]] <- list(fit=fit, pred=pred, actual=actual, resid=resid, mse=mse,
                        train_id=train_id, test_id=test_id)
  #})
  }
  
  # fit for the all data
  if (verbose) {
    message("Whole data fitting")
  }
  fit <- trlasso(X, y,
    beta_tilde = beta_tilde, family = family, lambda = lambda,
    alpha = alpha, penalty.factor.lasso = penalty.factor.lasso,
    penalty.factor.trlasso = penalty.factor.trlasso,
    penalty.factor.thresh = penalty.factor.thresh, eps = eps,
    verbose = verbose, ...
  )
  # pred <- predict.trlasso(fit, X)
  
  res <- list()
  res$lambda <- lambda
  pred_all <- Reduce(rbind, lapply(cv_fit, "[[", "pred"))
  actual_all <- Reduce("c", lapply(cv_fit, "[[", "actual"))

  resid_all <- cv.trlasso.resid_all(family, actual_all, pred_all)
  if (unit == "sample") {
    res$cve <- apply(resid_all, 2, mean)
    res$cvse <- apply(resid_all, 2, sd) / sqrt(nrow(resid_all))
  } else if (unit=="fold") {
    res$cve <- apply(sapply(cv_fit, "[[", "mse"), 1, mean)
    res$cvse <- apply(sapply(cv_fit, "[[", "mse"), 1, sd) / sqrt(nfolds)
  } else {
    stop(paste0(unit, "is not defined"))
  }
  res$cvup <- res$cve + res$cvse
  res$cvlo <- res$cve - res$cvse

  res$fit <- fit
  res$lambda.min.index <- which.min(res$cve)
  res$lambda.min <- lambda[res$lambda.min.index]
  res$lambda.1se.index <- suppressWarnings(
    max(Filter(function(i){
      i < res$lambda.min.index
    }, which(res$cve > res$cve[res$lambda.min.index] + res$cvse[res$lambda.min.index])))
  ) # -Inf if there are no lambda.1se.index
  res$lambda.1se <- res$lambda[res$lambda.1se.index]
  res$foldid <- foldid
  res$alpha <- alpha
  res$fit.preval <- pred_all
  res$actual <- actual_all
  
  class(res) <- "cv.trlasso"
  return(res)
}
