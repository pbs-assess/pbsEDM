% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{abs_mean_scaled_error}
\alias{abs_mean_scaled_error}
\title{Absolute Mean-Scaled Error}
\usage{
abs_mean_scaled_error(x, y)
}
\arguments{
\item{x}{[numeric()] [vector()] observed values}

\item{y}{[numeric()] [vector()] forecast values}
}
\value{
[numeric()] [vector()] absolute scaled errors
}
\description{
Absolute Mean-Scaled Error
}
\examples{
x <- seq(0, 10, 0.1)
y <- seq(10, 0, -0.1)
m <- abs_mean_scaled_error(x, y)
plot(m, type = "l")

}
