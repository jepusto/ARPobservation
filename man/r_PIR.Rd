% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_behavior_stream.R
\name{r_PIR}
\alias{r_PIR}
\title{Generates random partial interval recording behavior streams}
\usage{
r_PIR(
  n,
  mu,
  lambda,
  stream_length,
  F_event,
  F_interim,
  interval_length,
  rest_length = 0,
  summarize = FALSE,
  equilibrium = TRUE,
  p0 = 0,
  tuning = 2
)
}
\arguments{
\item{n}{number of behavior streams to generate}

\item{mu}{mean event duration}

\item{lambda}{mean interim time}

\item{stream_length}{length of behavior stream}

\item{F_event}{distribution of event durations. Must be of class \code{\link{eq_dist}}.}

\item{F_interim}{distribution of interim times. Must be of class \code{\link{eq_dist}}.}

\item{interval_length}{total interval length}

\item{rest_length}{length of any recording time in each interval}

\item{summarize}{logical value indicating whether the behavior streams should by summarized by taking their mean}

\item{equilibrium}{logical; if \code{TRUE}, then equilibrium initial conditions are used; 
if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
distributions are used for event durations and interim times.}

\item{p0}{Initial state probability. Only used if \code{equilibrium = FALSE}, in which case
default is zero (i.e., behavior stream always starts with an interim time).}

\item{tuning}{controls the size of the chunk of random event durations and interim times.
Adjusting this may be useful in order to speed computation time .}
}
\value{
If \code{summarize = FALSE}, a matrix with rows equal to \code{n} and a number of columns equal to the number intervals per session. If \code{summarize = TRUE} a vector of means of length \code{n}.
}
\description{
Random generation of behavior streams (based on an alternating
renewal process) of a specified length and with specified mean event 
durations, mean interim times, event distribution, and interim distribution,
which are then coded as partial interval recording data with given interval length
and rest length.
}
\details{
Generates behavior streams by repeatedly drawing random event durations and 
random interim times from the distributions as specified, until the sum of the durations and interim
times exceeds the requested stream length. Then applies a partial interval recording filter to the generated behavior streams.
}
\examples{

# An unsummarized set of PIR observations
r_PIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
       F_event = F_exp(), F_interim = F_exp(), 
       interval_length = 1, rest_length = 0)
      
# A summarized set of of PIR observations
r_PIR(n = 5, mu = 2, lambda = 4, stream_length = 20, 
       F_event = F_exp(), F_interim = F_exp(), 
       interval_length = 1, rest_length = 0,
       summarize = TRUE)
       
}
\author{
Daniel Swan <dswan@utexas.edu>
}
