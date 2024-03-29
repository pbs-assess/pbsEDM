% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbsEDM.R
\name{pbsEDM}
\alias{pbsEDM}
\title{Forecast via Empirical Dynamic Modelling}
\usage{
pbsEDM(
  N,
  lags,
  p = 1L,
  first_difference = FALSE,
  centre_and_scale = FALSE,
  exclusion_radius = "half",
  verbose = FALSE
)
}
\arguments{
\item{N}{A data frame with named columns for the response variable and
covariate time series.}

\item{lags}{A list of named integer vectors specifying the lags to use for
each time series in \code{N}.}

\item{p}{The integer forecast distance.}

\item{first_difference}{Logical. First-difference each time series?}

\item{centre_and_scale}{Logical. Centre and scale each time series?}

\item{exclusion_radius}{Number of points around ${\bf x}_t^*$ to exclude as
candidate nearest neighbours; either default of `half` as used for our manuscript
(see equation (6)), or a number to match
the `exclusionRadius` setting in `rEDM::Simplex()`. See `?pbsDist` for more details.}

\item{verbose}{Logical. Print progress?}
}
\value{
A list of class \code{pbsEDM} containing:

  \itemize{
\item \code{N} [matrix()] Response variable and unlagged covariates as
  columns

\item \code{N_observed} [vector()] Response variable time series

\item \code{N_forecast} [vector()] Forecast of response variable time series

\item \code{X} [matrix()] Unlagged and lagged state variables as columns

\item \code{X_observed} [vector()] Transformed response variable time series

\item \code{X_forecast} [vector()] Forecast of transformed response variable

\item \code{X_distance} [matrix()] Square distance \code{matrix} between
  pairs of points in state space (pairs of rows in \code{X})

\item \code{neighbour_distance} [matrix()] Distance by focal time (row) and
  rank (column)

\item \code{neighbour_index} [matrix()] Neighbour index by focal time (row)
  and distance rank (column)

\item \code{neighbour_value} [matrix()] Neighbour value by focal time (row)
  and distance rank (column)

\item \code{neighbour_weight} [matrix()] Neighbour weight by focal time (row)
  and distance rank (column)

\item \code{projected_index} [matrix()] Projected neighbour index by
  projected time (row) and neighbour distance rank (column)

\item \code{projected_value} [matrix()] Projected neighbour value by
  projected time (row) and neighbour distance rank (column)

\item \code{projected_weight} [matrix()] Projected neighbour weight by
  projected time (row) and neighbour distance rank (column)

\item \code{lags} [list()] A named list of integer vectors specifying the
  lags to use for each time  series in \code{N}

\item \code{p} [integer()] The forecast distance

\item \code{first_difference} [logical()] First difference each time series?

\item \code{centre_and_scale} [logical()] Centre and scale each time series?

\item \code{results} [data.frame()] A summary of forecast accuracy
}
}
\description{
Perform short-term nonlinear forecasting via Empirical Dynamic
  Modelling.
}
\details{
The name of the first element in \code{lags} must match the name of
  the response variable in \code{N}. Unlagged time series, including the
  response variable, must be specified by a zero in the corresponding named
  vector in \code{lags}. For example, given a \code{data.frame} with named
  columns \code{Predator}, \code{Prey} and \code{Temperature},
  \code{Predator} can be specified as the unlagged response variable by

  \code{lags = list(Predator = c(0, ...), ...)}.

  This places the unlagged time series of \code{Predator} abundance (or its
  optionally first-differenced and/or centred and scaled counterpart) along
  the first axis of the reconstructed state space. To predict \code{Predator}
  abundance from its first two lags, and from the unlagged and first lags of
  \code{Prey} and \code{Temperature}, \code{lags} can be specified as

  \code{lags = list(Predator = c(0:2), Prey = c(0:1), Temperature = c(0:1))}.

  This example generalizes to arbitrary (possibly non-consecutive) lags of
  arbitrarily many covariates (up to limitations of time series length).
}
\examples{
N <- matrix(rep(1:30, 5), ncol = 5)
colnames(N) <- c("A", "B", "C", "D", "E")
lags <- list(A = c(0, 1, 2), B = c(0, 1), C = c(0, 1, 2))
m1 <- pbsEDM(N, lags, verbose = TRUE)

N <- data.frame(x = simple_ts)
lags <- list(x = 0:1)
m2 <- pbsEDM(N, lags, verbose = TRUE)

N <- data.frame(x = simple_ts)
lags <- list(x = 0:1)
m3 <- pbsEDM(N, lags, first_difference = TRUE, verbose = TRUE)

}
\author{
Luke A. Rogers
}
