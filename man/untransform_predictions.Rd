% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/untransform_predictions.R
\name{untransform_predictions}
\alias{untransform_predictions}
\title{Untransform Vector of State Space Predictions}
\usage{
untransform_predictions(
  response_observed,
  response_s_predicted,
  max_lag,
  positive_response_only = TRUE
)
}
\arguments{
\item{response_observed}{[numeric()][vector()] original observed values}

\item{response_s_predicted}{[numeric()][vector()] state space predictions (from the
state space that used lagged and scaled variables), though just the
response column}

\item{max_lag}{maximum lag used}

\item{positive_response_only}{logical, if TRUE then set any negative
predictions of the response variable to be the minimum observed value}
}
\value{
[numeric()][vector()] a prediction vector of original `N_t` in
  absolute space, so
  `response_predicted`, but also with a value for `T+1`, which is what we are
  ultimately after.
}
\description{
So this is the reverse of `state_space_reconstruction_for_sve()`, but we only
  need to do it for the response variable, so just a single vector. Was using the
  notation of in Appendix of first manuscript, but only one vector so
  changing to `response_` to match the call from `single_view_embedding_for_sve()`.
}
