# library(testthat)
# source_withcheck <- function(path) {
#   results <- suppressWarnings(try(source(path), silent = TRUE)
#   )
#   if (class(results) == "try-error")
#   {
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
#   results <- suppressWarnings(try(sourceCpp(path), silent = TRUE)
#   )
#   if (class(results) == "try-error")
#   {
#     sourceCpp(paste0("../../", path))
#   }
# }
# sourceCpp_withcheck("./src/cda.cpp")

test_that("lambda sequence with gaussian", {
  n <- 30
  p <- 10
  
  set.seed(1234)
  
  x_test <-
    matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  beta_test <- c(2, 1, rep(0, length = p - 2))
  y_test <- x_test %*% beta_test + rnorm(n)
  
  lambda <- c(0.01, 0.1, 1)
  fit <- trlasso(
    X = x_test,
    y = y_test,
    lambda = lambda,
    family = "gaussian",
  )
  expect_equal(fit$lambda, lambda, tol = 1e-6)
  
  lambda <- c(
    1.55492403,
    0.93215246,
    0.55881072,
    0.33499823,
    0.20082616,
    0.12039212,
    0.07217318,
    0.04326668,
    0.02593770,
    0.01554924
  )
  fit <- trlasso(
    X = x_test,
    y = y_test,
    family = "gaussian",
    nlambda = 10
  )
  cv_fit <- cv.trlasso(
    X = x_test,
    y = y_test,
    family = "gaussian",
    nlambda = 10
  )
  expect_equal(fit$lambda, lambda, tol = 1e-6)
  expect_equal(cv_fit$lambda, lambda, tol = 1e-6)
})
