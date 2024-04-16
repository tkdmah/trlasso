# library(testthat)
# source_withcheck <- function(path) {
#   results <- suppressWarnings(try(source(path), silent = TRUE))
#   if (class(results) == "try-error") {
#     source(paste0("../../", path))
#   }
# }
# source_withcheck("./R/lasso.R")
# source_withcheck("./R/cross_validator.R")
# source_withcheck("./R/predictor.R")
# source_withcheck("./R/plotter.R")
# 
# library(Rcpp)
# sourceCpp_withcheck <- function(path) {
#   results <- suppressWarnings(try(sourceCpp(path), silent = TRUE))
#   if (class(results) == "try-error") {
#     sourceCpp(paste0("../../", path))
#   }
# }
# sourceCpp_withcheck("./src/cda.cpp")

test_that("cv.trlasso function with gaussian", {
  n <- 30
  p <- 10

  set.seed(1234)

  x_test <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  beta_test <- c(2, 1, rep(0, length = p - 2))
  y_test <- x_test %*% beta_test + rnorm(n)

  lambda <- c(1, 0.1, 0.01)
  cv_fit <- cv.trlasso(
    X = x_test,
    y = y_test,
    lambda = lambda,
    beta_tilde = matrix(rep(0, length = p * length(lambda)),
      nrow = p, ncol = length(lambda)
    ),
    family = "gaussian",
  )
  # cv_fit_glmnet <- cv.glmnet(
  #   x = x_test,
  #   y = y_test,
  #   lambda = lambda,
  #   family = "gaussian"
  # )
  # cv_fit_ncvreg <- cv.ncvreg(
  #   X = x_test,
  #   y = y_test,
  #   lambda = lambda,
  #   penalty = "lasso",
  #   family = "gaussian"
  # )
  
  beta_pred_old <- t(matrix(
    data = c(
      0.6302361, 1.704054, 1.95507221,
      0.0872549, 1.098350, 1.28368551,
      0.0000000, 0.000000, -0.21376068,
      0.0000000, 0.000000, 0.18810974,
      0.0000000, 0.000000, -0.17300145,
      0.0000000, 0.000000, -0.02985586,
      0.0000000, 0.000000, 0.02699009,
      0.0000000, 0.000000, -0.14231951,
      0.0000000, 0.000000,  0.22738957,
      0.0000000, 0.000000, 0.05371640
    ),
    nrow = 3,
    ncol = 10,
  ))

  # print(cv_fit$fit$beta)
  expect_equal(cv_fit$fit$beta, beta_pred_old, tol = 1e-6)
})
