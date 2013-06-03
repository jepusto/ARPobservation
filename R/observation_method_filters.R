

## continuous duration recording #### 

CDR_single <- function(b_stream, stream_length) {
  start_state <- b_stream$start_state
  switches <- length(b_stream$b_stream)
  if (switches > 1) {
    duration <- ((if (!start_state==(switches %% 2)) stream_length else 0) +
                   sum(b_stream$b_stream[seq(2 - start_state, switches, 2)]) -
                   sum(b_stream$b_stream[seq(1 - start_state, switches, 2)])) / stream_length
  } else if (switches == 1) {
    duration <- ((2 * start_state - 1) * b_stream$b_stream + 
                   (1 - start_state) * stream_length) / stream_length
  } else {
    duration <- start_state
  }
  return(duration)
}

#' @title Applies continuous duration recording to a behavior stream
#' 
#' @description
#' Calculates the proportion of session time during which behavior occurs.
#' 
#' @param BS object of class \code{behavior_stream}
#' 
#' @export 
#' 
#' @return Vector of proportions.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' continuous_duration_recording(BS)

continuous_duration_recording <- function(BS) 
  sapply(BS$b_streams, CDR_single, stream_length = BS$stream_length)





## momentary time recording ####

MTS_single <- function(b_stream, stream_length, n_intervals, summarize = TRUE) {
  obs_pos <- sapply(stream_length * (summarize:n_intervals) / n_intervals, 
                    findInterval, vec=b_stream$b_stream)
  MTS <- (obs_pos %% 2) != b_stream$start_state
  if (summarize) return(mean(MTS)) else return(MTS)
}

#' @title Applies momentary time recording to a behavior stream
#' 
#' @description Evaluates the presence or absence of the behavior at fixed moments in time.
#' 
#' @param BS object of class behavior_stream
#' @param n_intervals number of intervals into which the behavior stream is divided.
#' @param summarize logical value indicating whether vector of moments should be summarized by taking their mean.
#' 
#' @export 
#' 
#' @return If \code{summarize = FALSE}, a matrix with length \code{n_intervals + 1} and width equal to the number
#' of behavior streams in \code{BS}. If \code{summarize = TRUE}, a vector of proportions of length equal to the 
#' number of behavior streams in \code{BS}. Note that if \code{summarize = TRUE}, the initial state of the 
#' behavior stream is excluded when calculating the mean, so the proportion is based on \code{n_intervals} values.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' momentary_time_recording(BS, 30, FALSE)
#' momentary_time_recording(BS, 30)
#' colMeans(momentary_time_recording(BS, 30, FALSE)[-1,])

momentary_time_recording <- function(BS, n_intervals, summarize = TRUE) 
  sapply(BS$b_streams, MTS_single, 
         stream_length = BS$stream_length, 
         n_intervals = n_intervals, summarize = summarize)





## event counting ####

#' @title Applies event counting to a behavior stream
#' 
#' @description
#' Calculates the number of behaviors that begin during the observation session.
#' 
#' @param BS object of class \code{behavior_stream}
#' 
#' @export 
#' 
#' @return Vector of non-negative integers.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' event_counting(BS)

event_counting <- function(BS) 
  sapply(BS$b_streams, function(x) 
    floor((length(x$b_stream) + 1 - x$start_state) / 2))





## interval recording ####

IntRec_single <- function(b_stream, stream_length, 
                          n_intervals, rest_proportion = 0, 
                          partial = TRUE, summarize = TRUE) {
  int_start <- sapply(stream_length * (0:(n_intervals-1)) / n_intervals, findInterval, vec=b_stream$b_stream)
  int_end <- sapply(stream_length * ((1:n_intervals) - rest_proportion) / n_intervals, findInterval, vec=b_stream$b_stream)
  WIR <- (int_start == int_end) * ((int_start + 1 - partial) %% 2 == b_stream$start_state)
  if (summarize) return(ifelse(partial, 1 - mean(WIR), mean(WIR))) else
    if (partial) return(1 - WIR) else return(WIR)
}

#' @title Applies interval recording to a behavior stream
#' 
#' @description Divides the observation session into a specified number of intervals. For partial interval recording, 
#' each interval is scored according to whether the behavior is present at any point during the interval. For whole
#' interval recording, each interval is scored according to whether the behavior is present for the duration.
#' 
#' @param BS object of class behavior_stream
#' @param n_intervals number of intervals into which the behavior stream is divided.
#' @param rest_proportion proportion of each interval to exclude from observation. See details.
#' @param partial logical value indicating whether to use partial interval recording (\code{TRUE}) or 
#'        whole interval recording (\code{FALSE}).
#' @param summarize logical value indicating whether vector of moments should be summarized by taking their mean.
#' 
#' @details
#' Each behavior stream is divided into \code{n_intervals} intervals. The last \code{rest_proportion} of each interval is 
#' excluded from observation. For example, for a stream length of 100, \code{n_intervals = 20}, and 
#' \code{rest_proportion = 0.2}, the first interval runs from [0,4), the second interval runs from [5,9], etc. 
#' 
#' @export 
#' 
#' @return If \code{summarize = FALSE}, a matrix with length \code{n_intervals} and width equal to the number
#' of behavior streams in \code{BS}. If \code{summarize = TRUE}, a vector of proportions of length equal to the 
#' number of behavior streams in \code{BS}. 
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' interval_recording(BS, 30, partial = TRUE, summarize = FALSE)
#' interval_recording(BS, 30, partial = TRUE, summarize = TRUE)
#' colMeans(interval_recording(BS, 30, partial = TRUE, summarize = FALSE))
#' interval_recording(BS, 30, partial = FALSE)


interval_recording <- function(BS, n_intervals, rest_proportion = 0, partial = TRUE, summarize = TRUE)
  sapply(BS$b_streams, IntRec_single,
         stream_length = BS$stream_length, n_intervals = n_intervals, 
         rest_proportion = rest_proportion, partial = partial, summarize = summarize)





## Process behavior stream using multiple procedures ####

#' @title Applies multiple recording procedures to a behavior stream
#' 
#' @description This is a convenience function that allows multiple recording procedures to be applied
#' to a single behavior stream. Results are reported either per behavior stream or as summary statistics, averaged
#' over multiple behavior streams.
#' 
#' @param BS object of class behavior_stream
#' @param data_types list of recording procedures to apply to the behavior stream. See details.
#' @param n_intervals number of intervals into which the behavior stream is divided, for purposes of momentary
#'        time recording and interval recording procedures.
#' @param rest_proportion proportion of each interval to exclude from observation for interval recording. See 
#'        documentation for \code{\link{interval_recording}}.
#' @param n_aggregate number of observations over which to calculate summary statistics.
#' 
#' @details
#' The following recording procedures are currently implemented
#' \itemize{
#' \item \code{C} - continuous duration recording
#' \item \code{M} - momentary time recording
#' \item \code{E} - event counting
#' \item \code{P} - partial interval recording
#' \item \code{W} - whole interval recording
#' }
#' 
#' @export 
#' 
#' @return If \code{n_aggregate = 1}, a data frame with one column per procedure listed in \code{data_types} and length equal to the number
#' of behavior streams in \code{BS}. If \code{n_aggregate > 1}, a list containing two data frames: one with sample means
#' and one with sample variances, both taken across \code{n_aggregate} behavior streams.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 50, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' reported_observations(BS, n_intervals = 30)
#' reported_observations(BS, n_intervals = 30, n_aggregate = 5)

reported_observations <- function(BS, data_types = c("C","M","E","P","W"),
                                  n_intervals = 0, rest_proportion = 0, n_aggregate = 1) {
  n <- length(BS$b_streams)
  recorded <- data.frame(matrix(NA, n, length(data_types)))
  names(recorded) <- data_types
  
  if ("C" %in% data_types) 
    recorded$CDR <- continuous_duration_recording(BS)  
  if ("M" %in% data_types) 
    recorded$MTR <- momentary_time_recording(BS, n_intervals)
  if ("E" %in% data_types) 
    recorded$EC <- event_counting(BS)
  if ("P" %in% data_types) 
    recorded$PIR <- interval_recording(BS, n_intervals, rest_proportion, partial = TRUE)  
  if ("W" %in% data_types) 
    recorded$WIR <- interval_recording(BS, n_intervals, rest_proportion, partial = FALSE)  
  
  if (n_aggregate == 1) return(recorded) else {
    groups <- list(rep(1:(n / n_aggregate), each = n_aggregate))
    recorded_mean <- aggregate(recorded, groups, mean)
    recorded_var <- aggregate(recorded, groups, var)
    return(list(mean = recorded_mean[,-1], var = recorded_var[,-1]))
  }
}





## augmented interval recording ####

#' @title Applies augmented interval recording to a behavior stream
#' 
#' @description Divides the observation session into a specified number of intervals. 
#' Each interval is scored using partial interval recording, whole interval recording, and  
#' momentary time sampling (at the beginning of the following interval). The sum of the three
#' scores is then recorded.
#' 
#' @param BS object of class behavior_stream
#' @param n_intervals number of intervals into which the behavior stream is divided.
#' @param rest_proportion proportion of each interval to exclude from observation. See details.
#' 
#' @details
#' Each behavior stream is divided into \code{n_intervals} intervals. The last \code{rest_proportion} of each interval is 
#' excluded from observation. For example, for a stream length of 100, \code{n_intervals = 20}, and 
#' \code{rest_proportion = 0.2}, the first interval runs from [0,4), the second interval runs from [5,9], etc. 
#' 
#' @export 
#' 
#' @return A matrix with length \code{n_intervals} and width equal to the number
#' of behavior streams in \code{BS}.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' augmented_recording(BS, 30)

augmented_recording <- function(BS, n_intervals, rest_proportion = 0) {
  PIR <- interval_recording(BS, n_intervals, rest_proportion, partial = TRUE, summarize = FALSE)
  WIR <- interval_recording(BS, n_intervals, rest_proportion, partial = FALSE, summarize = FALSE)
  MTS <- momentary_time_recording(BS, n_intervals, summarize = FALSE)
  list(initial = MTS[1,], AIR = PIR + WIR + MTS[-1,])
}





## intermittent transition recording ####

#' @title Applies intermittent transition recording to a behavior stream
#' 
#' @description Divides the observation session into a specified number of intervals. 
#' Each interval is scored using partial interval recording, whole interval recording, and  
#' momentary time sampling (at the beginning of the following interval). The sum of the three
#' scores is then recorded.
#' 
#' @param BS object of class behavior_stream
#' @param n_intervals number of intervals into which the behavior stream is divided.
#' @param rest_proportion proportion of each interval to exclude from observation. See details.
#' 
#' @details
#' Each behavior stream is divided into \code{n_intervals} intervals. The last \code{rest_proportion} of each interval is 
#' excluded from observation. For example, for a stream length of 100, \code{n_intervals = 20}, and 
#' \code{rest_proportion = 0.2}, the first interval runs from [0,4), the second interval runs from [5,9], etc. 
#' 
#' @export 
#' 
#' @return A matrix with length \code{n_intervals} and width equal to the number
#' of behavior streams in \code{BS}.
#' 
#' @examples
#' BS <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                        F_mu = F_exp(), F_lambda = F_exp(), stream_length = 100)
#' augmented_recording(BS, 30)

intermittent_transition_recording <- function(BS, n_intervals, rest_proportion = 0) {
  PIR <- interval_recording(BS, n_intervals, rest_proportion, partial = TRUE, summarize = FALSE)
  WIR <- interval_recording(BS, n_intervals, rest_proportion, partial = FALSE, summarize = FALSE)
  MTS <- momentary_time_recording(BS, n_intervals, summarize = FALSE)
  list(initial = MTS[1,], AIR = PIR + WIR + MTS[-1,])
}
