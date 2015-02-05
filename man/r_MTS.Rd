% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/simulate_behavior_stream.R
\name{r_MTS}
\alias{r_MTS}
\title{Generates random momentary time sampling behavior streams}
\usage{
r_MTS(n, mu, lambda, stream_length, F_event, F_interim, interval_length,
  equilibrium = TRUE, p0 = 0, tuning = 2)
}
\arguments{
\item{n}{number of behavior streams to generate}

\item{mu}{mean event duration}

\item{lambda}{mean interim time}

\item{stream_length}{length of behavior stream}

\item{F_event}{distribution of event durations. Must be of class \code{\link{eq_dist}}.}

\item{F_interim}{distribution of interim times. Must be of class \code{\link{eq_dist}}.}

\item{interval_length}{length of time between moments}

\item{equilibrium}{logical; if \code{TRUE}, then equilibrium initial conditions are used;
if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating
distributions are used for event durations and interim times.}

\item{p0}{vector of initial state probabilities. Only used if \code{equilibrium = FALSE}, in which case
default is zero (i.e., behavior stream always starts with an interim time).}

\item{tuning}{controls the size of the chunk of random event durations and interim times.
Adjusting this may be useful in order to speed computation time .}
}
\value{
A matrix with rows equal to \code{n} and length equal to \code{(stream_length/interval_length) + 1}
}
\description{
Random generation of behavior streams (based on an alternating
renewal process) of a specified length and with specified mean event
durations, mean interim times, event distribution, and interim distribution,
which are then coded as momentary time sampling data with given interval length
between moments.
}
\details{
Generates behavior streams by repeatedly drawing random event durations and
random interim times from the distributions as specified, until the sum of the durations and interim
times exceeds the requested stream length. Then applies a momentary time sampling filter to the generated behavior streams.
}
