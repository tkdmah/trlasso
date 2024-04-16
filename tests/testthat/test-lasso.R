# library(testthat)
# source_withcheck <- function(path) {
#   results <- suppressWarnings(try(source(path), silent = TRUE))
#   if (class(results) == "try-error") {
#     source(paste0("../../", path))
#   }
# }
# source_withcheck("./R/lasso.R")
# 
# library(Rcpp)
# sourceCpp_withcheck <- function(path) {
#   results <- suppressWarnings(try(sourceCpp(path), silent = TRUE))
#   if (class(results) == "try-error") {
#     sourceCpp(paste0("../../", path))
#   }
# }
# sourceCpp_withcheck("./src/cda.cpp")

test_that("trlasso function with gaussian", {
  n <- 10
  p <- 10

  set.seed(1234)

  x_test <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  beta_test <- c(2, 1, rep(0, length = p - 2))
  y_test <- x_test %*% beta_test + rnorm(n)

  lambda <- c(1, 0.1, 0.01)
  fit <- trlasso(
    X = x_test,
    y = y_test,
    lambda = lambda,
    beta_tilde = matrix(rep(0, length = p * length(lambda)),
      nrow = p, ncol = length(lambda)
    ),
    family = "gaussian",
  )
  # fit_ncvreg <- ncvreg(
  #   X = x_test,
  #   y = y_test,
  #   lambda = lambda,
  #   penalty = "lasso",
  #   family = "gaussian"
  # )
  # fit_glmnet <- glmnet(
  #   x = x_test,
  #   y = y_test,
  #   lambda = lambda,
  #   family = "gaussian"
  # )
  
  beta_pred_old <- t(matrix(
    data = c(
      0.6929767, 1.72589715, 1.88870083,
      0.0000000,  0.65100672, 0.72893186,
      -0.1675009, -0.43551455, -0.66260901,
      0.0000000, 0.01946874, 0.00000000,
      0.0000000, 0.11012248, 0.52985515,
      0.0000000, 0.00000000, 0.19824745,
      0.0000000, 0.00000000, -0.03436448,
      0.0000000, 0.00000000, 0.00000000,
      0.0000000, 0.00000000, -0.19380120,
      0.0000000, 0.13596144, 0.15577425
    ),
    nrow = 3,
    ncol = 10,
  ))

  # print(fit)
  expect_equal(fit$beta, beta_pred_old, tol = 1e-6)
})
