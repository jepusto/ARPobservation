\name{reported_observations}
\alias{reported_observations}
\title{Applies multiple recording procedures to a behavior stream}
\usage{
  reported_observations(BS,
    data_types = c("C", "M", "E", "P", "W"),
    n_intervals = 0, rest_proportion = 0, n_aggregate = 1)
}
\arguments{
  \item{BS}{object of class behavior_stream}

  \item{data_types}{list of recording procedures to apply
  to the behavior stream. See details.}

  \item{n_intervals}{number of intervals into which the
  behavior stream is divided, for purposes of momentary
  time recording and interval recording procedures.}

  \item{rest_proportion}{proportion of each interval to
  exclude from observation for interval recording. See
  documentation for \code{\link{interval_recording}}.}

  \item{n_aggregate}{number of observations over which to
  calculate summary statistics.}
}
\value{
  If \code{n_aggregate = 1}, a data frame with one column
  per procedure listed in \code{data_types} and length
  equal to the number of behavior streams in \code{BS}. If
  \code{n_aggregate > 1}, a list containing two data
  frames: one with sample means and one with sample
  variances, both taken across \code{n_aggregate} behavior
  streams.
}
\description{
  This is a convenience function that allows multiple
  recording procedures to be applied to a single behavior
  stream. Results are reported either per behavior stream
  or as summary statistics, averaged over multiple behavior
  streams.
}
\details{
  The following recording procedures are currently
  implemented \itemize{ \item \code{C} - continuous
  duration recording \item \code{M} - momentary time
  recording \item \code{E} - event counting \item \code{P}
  - partial interval recording \item \code{W} - whole
  interval recording }
}
\examples{
BS <- r_behavior_stream(n = 50, mu = 3, lambda = 10,
                       F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
reported_observations(BS, n_intervals = 30)
reported_observations(BS, n_intervals = 30, n_aggregate = 5)
}

