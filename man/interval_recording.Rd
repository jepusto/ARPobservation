\name{interval_recording}
\alias{interval_recording}
\title{Applies interval recording to a behavior stream}
\usage{
  interval_recording(BS, n_intervals, rest_proportion = 0,
    partial = TRUE, summarize = TRUE)
}
\arguments{
  \item{BS}{object of class behavior_stream}

  \item{n_intervals}{number of intervals into which the
  behavior stream is divided.}

  \item{rest_proportion}{proportion of each interval to
  exclude from observation. See details.}

  \item{partial}{logical value indicating whether to use
  partial interval recording (\code{TRUE}) or whole
  interval recording (\code{FALSE}).}

  \item{summarize}{logical value indicating whether vector
  of moments should be summarized by taking their mean.}
}
\value{
  If \code{summarize = FALSE}, a matrix with length
  \code{n_intervals} and width equal to the number of
  behavior streams in \code{BS}. If \code{summarize =
  TRUE}, a vector of proportions of length equal to the
  number of behavior streams in \code{BS}.
}
\description{
  Divides the observation session into a specified number
  of intervals. For partial interval recording, each
  interval is scored according to whether the behavior is
  present at any point during the interval. For whole
  interval recording, each interval is scored according to
  whether the behavior is present for the duration.
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
interval_recording(BS, 30, partial = TRUE, summarize = FALSE)
interval_recording(BS, 30, partial = TRUE, summarize = TRUE)
colMeans(interval_recording(BS, 30, partial = TRUE, summarize = FALSE))
interval_recording(BS, 30, partial = FALSE)
}

