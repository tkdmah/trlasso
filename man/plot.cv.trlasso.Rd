% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotter.R
\name{plot.cv.trlasso}
\alias{plot.cv.trlasso}
\title{Plot a cross validation error path}
\usage{
\method{plot}{cv.trlasso}(x, ...)
}
\arguments{
\item{x}{cross validated trlasso model}

\item{...}{parameters of}
}
\description{
Plot a cross validation error path
}
\examples{
X <- as.matrix(iris[, 1:3])
y <- iris[, 4]
cv_fit <- cv.trlasso(X, y, nlambda=50, lambda.min.ratio=1e-2)
plot(cv_fit)
plot(cv_fit$fit)

}
