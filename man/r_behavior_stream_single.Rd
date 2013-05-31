\name{r_behavior_stream_single}
\alias{r_behavior_stream_single}
\title{Generate a single random behavior stream}
\usage{
  r_behavior_stream_single(mu, lambda, F_mu, F_lambda,
    stream_length, tuning = 2)
}
\arguments{
  \item{mu}{mean event duration}

  \item{lambda}{mean interim time}

  \item{F_mu}{distribution of event durations. Must be of
  class \code{\link{eq_dist}}.}

  \item{F_lambda}{distribution of interim times. Must be of
  class \code{\link{eq_dist}}.}

  \item{stream_length}{length of behavior stream}

  \item{tuning}{controls the size of the chunk of random
  event durations and interim times. Adjusting this may be
  useful in order to speed computation time .}
}
\value{
  Returns a list containing the starting state and the
  vector of transition times between states.
}
\description{
  Random generation of a single behavior stream (based on
  an alternating renewal process) of a specified length and
  with specified mean event duration, mean interim time,
  event distribution, and interim distribution.
}
\details{
  Generates a single random behavior stream by repeatedly
  drawing random event durations and random interim times
  from the distributions as specified, until the sum of the
  durations and interim times exceeds the requested stream
  length.
}

