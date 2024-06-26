% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotter.R
\name{plot.trlasso}
\alias{plot.trlasso}
\title{Plot a solution path}
\usage{
\method{plot}{trlasso}(x, ...)
}
\arguments{
\item{x}{trlasso model}

\item{...}{parameters of matlines function}
}
\description{
Plot a solution path
}
\examples{
X <- as.matrix(iris[, 1:3])
y <- iris[, 4]
fit <- trlasso(X, y, nlambda=50, lambda.min.ratio=1e-2)
plot(fit)

}
