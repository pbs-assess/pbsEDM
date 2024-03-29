% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_time_series}
\alias{plot_time_series}
\title{Plot the observed time series as either `N_t` or `Y_t`}
\usage{
plot_time_series(
  values,
  X.or.N,
  par.mar.ts = c(3, 3, 1, 1),
  max_time = NULL,
  t.axis.range = NULL,
  last.time.to.plot = NULL,
  late.num = 3,
  late.col = "red",
  early.col = "black",
  early.col.lines = "lightgrey",
  start = 1,
  pt.type = "p",
  par.mgp = c(1.5, 0.5, 0),
  add.legend = TRUE,
  label = NULL,
  label.cex = 0.7
)
}
\arguments{
\item{values}{vector of values to be plotted}

\item{X.or.N}{"N" if raw non-differenced data, "X" for differenced data}

\item{par.mar.ts}{`par(mar)` values}

\item{max_time}{maximum time value for the time axis}

\item{t.axis.range}{range of time axis}

\item{last.time.to.plot}{last time value of N[t] to use when plotting, so
final Y[t] used will be Y[t-1] (since Y[t] uses N[t+1])}

\item{late.num}{final number of `N[t]` time steps to plot in a different colour}

\item{late.col}{colour in which to plot final `late.num` time steps}

\item{early.col}{colour in which to plot earlier time step points}

\item{early.col.lines}{colour in which to plot earlier time step points}

\item{start}{first time step (must be 1)}

\item{pt.type}{`type` value for `points()`}

\item{par.mgp}{`par("mgp")` values}

\item{add.legend}{logical, whether to add legend to the `N[t]` time series}

\item{label}{label to annotate plot, such as `(a)` etc. for six-panel figure}

\item{label.cex}{size of label annotation}
}
\description{
First value must be `t=1`. For non-differenced values `N_t`, shows the time
series with the final values in a different colour, and a title showing the
final time step. For first-differenced values `Y_t`, shows the time series
of those, plus a one-dimensional phase plot.
}
\examples{
\donttest{
  plot_time_series(NY_lags_example$N_t, X.or.N = "N")
  plot_time_series(NY_lags_example$Y_t, X.or.N = "X")
}
}
\author{
Andrew Edwards
}
