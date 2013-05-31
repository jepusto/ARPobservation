
## generate single behavior stream ####

#' @title Generate a single random behavior stream
#' 
#' @description
#' Random generation of a single behavior stream (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' duration, mean interim time, event distribution, and interim distribution.
#' 
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param F_mu distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_lambda distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates a single random behavior stream by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. 
#' 
#' @return Returns a list containing the starting state and the vector of transition times between states.

r_behavior_stream_single <- function(mu, lambda, F_mu, F_lambda, stream_length, tuning = 2) {
  # initial condition
  start_state <- rbinom(1, 1, mu / (mu + lambda))
  
  # draw initial time
  if (start_state) {
    b_stream <- F_mu$r_eq(1, mu)
  } else {
    b_stream <- F_lambda$r_eq(1, lambda)
  }
  
  cum_length <- b_stream 
  cum_size <- 1
  
  # add event durations and interim times until total length exceeds stream length
  while (cum_length < stream_length) {
    
    # generate random event durations and interim times
    extend_size <- ceiling(tuning * (stream_length - cum_length) / (mu + lambda))
    event_times <- F_mu$r_gen(n=extend_size, mean = mu)
    interim_times <- F_lambda$r_gen(n=extend_size, mean = lambda)
    
    # lengthen behavior stream vector
    b_stream <- append(b_stream, cum_length + cumsum(
      if (start_state) c(rbind(interim_times, event_times)) 
      else c(rbind(event_times, interim_times))))
    
    # update totals
    cum_size <- cum_size + 2 * extend_size
    cum_length <- b_stream[cum_size]
  }
  
  list(start_state=start_state, b_stream = b_stream[b_stream < stream_length])
}



## generate behavior streams ####

#' @title Generates random behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution.
#' 
#' @param n number of behavior streams to generate
#' @param mu vector of mean event durations
#' @param lambda vector of mean interim time
#' @inheritParams r_behavior_stream_single
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. The vectors \code{mu} and \code{lambda} are recycled to length \code{n}.
#' \code{r_behavior_stream} works by calling \code{\link{r_behavior_stream_single}} repeatedly.
#' @export
#' 
#' @return An object of class \code{behavior_stream} containing two elements.
#' 
#' @examples
#' r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                   F_mu = F_exp(), F_lambda = F_exp(), 
#'                   stream_length = 100)

r_behavior_stream <- function(n, mu, lambda, F_mu, F_lambda, stream_length, tuning = 2) {
  
  mu_vec <- rep(mu, length.out = n)
  lambda_vec <- rep(lambda, length.out = n)
  
  BS <- list(stream_length = stream_length,
       b_streams = mapply(r_behavior_stream_single, 
                          mu = mu_vec,
                          lambda = lambda_vec,
                          MoreArgs = list(F_mu, F_lambda, stream_length, tuning),
                          SIMPLIFY = FALSE))
  class(BS) <- "behavior_stream"
  return(BS)
}


