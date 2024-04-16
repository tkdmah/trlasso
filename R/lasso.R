#' Fit a model using a design matrix
#'
#' @param X matrix of explanatory variables
#' @param y vector of objective variable
#' @param beta_tilde vector of initial parameters
#' @param family family of regression: "gaussian" is only supported
#' @param lambda.min.ratio ratio of max lambda and min lambda (ignored if lambda is specified)
#' @param nlambda the number of lambda (ignored if lambda is specified)
#' @param lambda lambda sequence
#' @param alpha alpha
#' @param eps convergence threshold (default 1e-4)
#' @param verbose whether output verbose warnings and messages (default FALSE)
#' @param standardize standardize or not (experimental, not safe)
#' @param penalty.factor.lasso factors applied to each coefficient in lasso penalty
#' @param penalty.factor.trlasso factors applied to each coefficient in trlasso penalty
#' @param penalty.factor.thresh lower limit for penalty factors (default 1e-2)
#' @param ... parameters for optimization
#'
#' @return lasso model
#' \item{beta}{coefficients}
#' \item{beta_standard}{standardized coefficients}
#' \item{a0}{intercepts}
#' \item{lambda}{regularization parameters}
#' \item{family}{family}
#'
#' @examples
#' X <- as.matrix(iris[, 1:3])
#' y <- iris[, 4]
#' fit <- trlasso(X, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(fit)
#'
#' @export
trlasso <- function(X, y, beta_tilde = NULL, family="gaussian", 
                    lambda.min.ratio =  1e-2,
                    nlambda = 100,
                    lambda = NULL,
                    alpha = 1,
                    eps = 1e-4,
                    verbose = FALSE,
                    standardize = TRUE,
                    penalty.factor.lasso = NULL,
                    penalty.factor.trlasso = NULL,
                    penalty.factor.thresh = 1e-2,
                    ...) {
  
  n <- nrow(X)

  t1 <- Sys.time()
  
  # constract standardized X and y
  if (verbose) {
    message("observed preprocessing")
  }
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
  y_sd <- sqrt(mean(y_tilde^2, na.rm=TRUE))
  
  # construct covariance matrices
  # Gamma <- (n-1) / n * cov(X_tilde, use=use) #sample covariance
  Gamma <- t(X_tilde) %*% X_tilde / n
  gamma <- t(X_tilde) %*% y_tilde / n #sample covariance
    
  t2 <- Sys.time()
  if (verbose) {
    message(paste0("Covariance matrix estimation takes ", format(t2-t1, digits = 3), " ", units(t2-t1)))
  }
  
  # initialize penalty
  if (is.null(penalty.factor.lasso)) {
    penalty.factor.lasso <- rep(1, length=ncol(X))
  } else {
    penalty.factor.lasso[penalty.factor.lasso < penalty.factor.thresh] <- penalty.factor.thresh
  }
  if  (is.null(penalty.factor.trlasso)) {
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
  
  # fit
  if (family == "gaussian") {
    if (is.null(beta_tilde)) {
      beta_tilde_standard <- rep(0, length = ncol(X))
    } else {
      beta_tilde_standard <- beta_tilde * X_sd
    }
    eps <- y_sd * eps
    if (verbose) message("start cov_lasso")
    fit <- cov_lasso(Gamma, gamma, as.numeric(beta_tilde_standard),
                     lambda = lambda, alpha = alpha, eps = eps,
                     penalty.factor.lasso = penalty.factor.lasso, penalty.factor.trlasso = penalty.factor.trlasso, ...
    )
    if (is.null(dim(fit$beta_standard))) {
      fit$beta <- fit$beta_standard / X_sd
      if (is.null(names(fit$beta))) {
        if (is.null(colnames(X))) {
          names(fit$beta) <- paste0("V", 1:length(fit$beta))
        } else {
          names(fit$beta) <- colnames(X)
        }
      }
    } else {
      fit$beta <- sweep(fit$beta_standard, 1, X_sd, FUN = "/")
      if (is.null(rownames(fit$beta))) {
        if (is.null(colnames(X))) {
          # rownames(fit$beta) <- paste0("V",1:length(fit$beta))
        } else {
          rownames(fit$beta) <- colnames(X)
        }
      }
    }
    fit$a0 <- mean(y) - t(X_mean) %*% fit$beta
  } else {
    stop("Specified family is not supported.")
  }

  t3 <- Sys.time()
  if (verbose) {
    message(paste0("Regression takes ", format(t3-t2, digits = 3), " ", units(t3-t2)))
  }
  
  fit$family <- family
  fit$X_mean <- X_mean
  fit$X_sd <- X_sd
  if (exists("Gamma")) {
    fit$Gamma <- Gamma
  } else {
    fit$Gamma <- 0
  }
  if (exists("gamma")) {
    fit$gamma <- gamma
  } else {
    fit$gamma <- 0
  }
  fit$alpha <- alpha
  class(fit) <- "trlasso"
  
  return(fit)
}


cov_lasso <- function(Gamma, gamma, beta_tilde = NULL,
                      lambda.min.ratio = 0.0001,
                      nlambda = 100,
                      lambda = NULL,
                      alpha = 1,
                      maxit = 1e+4,
                      eps = 1e-4,
                      init.beta = NULL,
                      strong = TRUE,
                      penalty.factor.lasso = NULL,
                      penalty.factor.trlasso = NULL,
                      penalty.factor.thresh = 1e-2) {
  
  # initialize lambda
  if (is.null(lambda)) {
    lambda_max <- max(gamma)
    lambda <- lambda_max * exp(seq(from=0,
                                   to=log(lambda.min.ratio),
                                   by=log(lambda.min.ratio)/(nlambda-1)))
  } else {
    nlambda <- length(lambda)
  }
  
  # initialize beta_tilde
  if (is.null(beta_tilde)) {
    beta_tilde <- rep(0, length(gamma))
  }
  
  # initialize warm
  init.beta <- matrix(rep(0, length=length(gamma)*length(lambda)),
                      nrow=length(gamma), ncol=length(lambda))
  
  # message("start covCdaC")
  beta <- covCdaC(Gamma, as.numeric(gamma), 
                  lambda, init.beta, as.numeric(beta_tilde), 
                  penalty.factor.lasso, penalty.factor.trlasso, alpha, 
                  maxit, eps, strong)
  
  # set result
  res <- list()
  res$beta_standard <- beta
  res$lambda <- lambda
  
  return(res)
}
