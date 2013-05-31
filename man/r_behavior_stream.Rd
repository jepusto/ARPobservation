\name{r_behavior_stream}
\alias{r_behavior_stream}
\title{Generates random behavior streams}
\usage{
  r_behavior_stream(n, mu, lambda, F_mu, F_lambda,
    stream_length, tuning = 2)
}
\arguments{
  \item{n}{number of behavior streams to generate}

  \item{mu}{vector of mean event durations}

  \item{lambda}{vector of mean interim time}

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
  An object of class \code{behavior_stream} containing two
  elements.
}
\description{
  Random generation of behavior streams (based on an
  alternating renewal process) of a specified length and
  with specified mean event durations, mean interim times,
  event distribution, and interim distribution.
}
\details{
  Generates behavior streams by repeatedly drawing random
  event durations and random interim times from the
  distributions as specified, until the sum of the
  durations and interim times exceeds the requested stream
  length. The vectors \code{mu} and \code{lambda} are
  recycled to length \code{n}. \code{r_behavior_stream}
  works by calling \code{\link{r_behavior_stream_single}}
  repeatedly.
}
\examples{
r_behavior_stream(n = 5, mu = 3, lambda = 10,
                  F_mu = F_exp(), F_lambda = F_exp(),
                  stream_length = 100)
}

