% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{rmse}
\alias{rmse}
\title{Root Mean Squared Error}
\usage{
rmse(x, y)
}
\arguments{
\item{x}{[numeric()][vector()] observed values}

\item{y}{[numeric()][vector()] forecast values}
}
\value{
[numeric()] root mean squared error
}
\description{
Root Mean Squared Error
}
\examples{
rmse(1:10, 10:1)

}
