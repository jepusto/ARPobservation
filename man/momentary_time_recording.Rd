\name{momentary_time_recording}
\alias{momentary_time_recording}
\title{Applies momentary time recording to a behavior stream}
\usage{
  momentary_time_recording(BS, n_intervals,
    summarize = TRUE)
}
\arguments{
  \item{BS}{object of class behavior_stream}

  \item{n_intervals}{number of intervals into which the
  behavior stream is divided.}

  \item{summarize}{logical value indicating whether vector
  of moments should be summarized by taking their mean.}
}
\value{
  If \code{summarize = FALSE}, a matrix with length
  \code{n_intervals + 1} and width equal to the number of
  behavior streams in \code{BS}. If \code{summarize =
  TRUE}, a vector of proportions of length equal to the
  number of behavior streams in \code{BS}. Note that if
  \code{summarize = TRUE}, the initial state of the
  behavior stream is excluded when calculating the mean, so
  the proportion is based on \code{n_intervals} values.
}
\description{
  Evaluates the presence or absence of the behavior at
  fixed moments in time.
}
\examples{
BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10,
                       F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
momentary_time_recording(BS, 30, FALSE)
momentary_time_recording(BS, 30)
colMeans(momentary_time_recording(BS, 30, FALSE)[-1,])
}

