\name{intermittent_transition_recording}
\alias{intermittent_transition_recording}
\title{Applies intermittent transition recording to a behavior stream}
\usage{
  intermittent_transition_recording(BS, n_intervals,
    rest_proportion = 0)
}
\arguments{
  \item{BS}{object of class behavior_stream}

  \item{n_intervals}{number of intervals into which the
  behavior stream is divided.}

  \item{rest_proportion}{proportion of each interval to
  exclude from observation. See details.}
}
\value{
  A matrix with length \code{n_intervals} and width equal
  to the number of behavior streams in \code{BS}.
}
\description{
  Divides the observation session into a specified number
  of intervals. Each interval is scored using partial
  interval recording, whole interval recording, and
  momentary time sampling (at the beginning of the
  following interval). The sum of the three scores is then
  recorded.
}
\details{
  Each behavior stream is divided into \code{n_intervals}
  intervals. The last \code{rest_proportion} of each
  interval is excluded from observation. For example, for a
  stream length of 100, \code{n_intervals = 20}, and
  \code{rest_proportion = 0.2}, the first interval runs
  from [0,4), the second interval runs from [5,9], etc.
}
\examples{
BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10,
                       F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
augmented_recording(BS, 30)
}

