% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{cv.trlasso.resid_all}
\alias{cv.trlasso.resid_all}
\title{Compute errors of predicted target variable.}
\usage{
cv.trlasso.resid_all(family, actual_all, pred_all, direct_prediction = FALSE)
}
\arguments{
\item{family}{family of regression: "gaussian"}

\item{actual_all}{matrix of actual target variable}

\item{pred_all}{matrix of predicted target variable}

\item{direct_prediction}{either corrected cross validation is used or not}
}
\description{
Compute errors of predicted target variable.
}
