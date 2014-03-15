
## generate single behavior stream ####

r_behavior_stream_single <- function(mu, lambda, F_mu, F_lambda, stream_length,
                                     equilibrium, p0, tuning) {
  # initial condition
  start_state <- rbinom(1, 1, p0)
  
  # draw initial time
  if (equilibrium) {
    if (start_state) {
      b_stream <- F_mu$r_eq(1, mu)
    } else {
      b_stream <- F_lambda$r_eq(1, lambda)
    } 
  } else {
    if (start_state) {
      b_stream <- F_mu$r_gen(1, mu)
    } else {
      b_stream <- F_lambda$r_gen(1, lambda)
    }    
    
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
#' @param F_mu distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_lambda distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param equilibrium logical; if \code{TRUE}, then equilibrium initial conditions are used; 
#' if \code{FALSE}, then \code{p0} is used to determine initial state and normal generating 
#' distributions are used for event durations and interim times.
#' @param p0 vector of initial state probabilities. Only used if \code{equilibrium = FALSE}, in which case
#' default is zero (i.e., behavior stream always starts with an interim time).
#' @param tuning controls the size of the chunk of random event durations and interim times.
#' Adjusting this may be useful in order to speed computation time .
#' 
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length. The vectors \code{mu}, \code{lambda}, and \code{p0} are 
#' recycled to length \code{n}.
#' @export
#' 
#' @return An object of class \code{behavior_stream} containing two elements.
#' 
#' @examples
#' # default equilibrium initial conditions
#' r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                   F_mu = F_exp(), F_lambda = F_exp(), 
#'                   stream_length = 100)
#'                   
#' # non-equilibrium initial conditions
#' r_behavior_stream(n = 5, mu = 3, lambda = 10,
#'                   F_mu = F_gam(3), F_lambda = F_gam(3),
#'                   stream_length = 100, 
#'                   equilibrium = FALSE, p0 = 0.5)

r_behavior_stream <- function(n, mu, lambda, F_mu, F_lambda, stream_length, 
                              equilibrium = TRUE, p0 = 0, tuning = 2) {
  
  mu_vec <- rep(mu, length.out = n)
  lambda_vec <- rep(lambda, length.out = n)
  p0_vec <- if (equilibrium) mu_vec / (mu_vec + lambda_vec) else rep(p0, length.out = n)
  
  BS <- list(stream_length = stream_length,
       b_streams = mapply(r_behavior_stream_single, 
                          mu = mu_vec,
                          lambda = lambda_vec,
                          p0 = p0_vec,
                          MoreArgs = list(F_mu = F_mu, F_lambda = F_lambda, 
                                      stream_length = stream_length, 
                                      equilibrium = equilibrium,
                                      tuning = tuning),
                          SIMPLIFY = FALSE))
  class(BS) <- "behavior_stream"
  return(BS)
}


