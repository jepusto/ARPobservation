% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effect_size_estimators.R
\name{PIR_MOM}
\alias{PIR_MOM}
\title{Moment estimator for prevalence and incidence, with bootstrap confidence intervals}
\usage{
PIR_MOM(
  PIR,
  phase,
  base_level,
  intervals,
  interval_length,
  rest_length = 0,
  Bootstraps = 2000,
  conf_level = 0.95,
  exponentiate = FALSE,
  seed = NULL
)
}
\arguments{
\item{PIR}{vector of PIR measurements}

\item{phase}{factor or vector indicating levels of the PIR measurements.}

\item{base_level}{a character string or value indicating the name of the baseline level.}

\item{intervals}{the number of intervals in the sample of observations}

\item{interval_length}{the total length of each interval}

\item{rest_length}{length of the portion of the interval devoted to recording. Default is \code{0}}

\item{Bootstraps}{desired number of bootstrap replicates. Default is \code{2000}}

\item{conf_level}{Desired coverage rate of the calculated confidence interval. Default is \code{.95}.}

\item{exponentiate}{a logical indicating whether the row corresponding to the ratio of treatment to baseline should be exponentiated, with the default as \code{FALSE}.}

\item{seed}{seed value set in order to make bootstrap results reproducible. Default is \code{null}}
}
\value{
A dataframe with six columns and three rows corresponding to baseline, treatment,
and the log ratio or ratio (depending upon the value of \code{exponentiate}) of treatment to baseline
}
\description{
Estimates prevalence and incidence for two samples,
along with the ratios of each parameter, assuming that the behavior follows
an `. Also provides bootstrap confidence intervals.
}
\details{
The moment estimators are based on the assumption that the
underlying behavior stream follows an Alternating Poisson Process, in which both the
event durations and interim times are exponentially distributed.

The \code{PIR} vector can be in any order corresponding to the factor or vector
\code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B",
"base" and "treat", or "0" and "1". If there are more than two levels in \code{phase}
this function will not work. A value for \code{base_level} must be specified - if it is a
character string it is case sensitive.

\code{intervals}, \code{interval_length}, and \code{rest_length} are all single values that
are assumed to be held constant across both samples and all observation sessions.
If vectors of values are provided for these variables, it is assumed that the first value
in each vector is constant across all observations.

\code{interval_length} This is the total length of each individual interval.
Sometimes a portion of the interval is set aside for recording purposes, in which case
\code{rest_length} should be set to the length of time devoted to recording.
The default assumption is that there is no recording time.
The length of time devoted to active observation is calculated as \code{interval_length - rest_length}.

At the default setting of \code{bootstraps = 2000}, PIR_MOM takes just under six seconds to run on an Intel Core i5-2410M processor.
}
\examples{

# Estimate prevalence and incidence ratios for Carl from the Moes dataset
data(Moes)
with(subset(Moes, Case == "Carl"),
  PIR_MOM(PIR = outcome,
          phase = Phase,
          intervals = intervals,
          interval_length = (active_length + rest_length),
          rest_length = rest_length,
          base_level = "No Choice",
          Bootstraps = 200,
          seed = 149568373))

}
\author{
Daniel Swan <dswan@utexas.edu>
}
