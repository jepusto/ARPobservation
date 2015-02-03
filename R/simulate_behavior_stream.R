
## generate single behavior stream ####

r_behavior_stream_single <- function(mu, lambda, F_event, F_interim, stream_length,
                                     equilibrium, p0, tuning) {
  # initial condition
  start_state <- rbinom(1, 1, p0)
  
  # draw initial time
  if (equilibrium) {
    if (start_state) {
      b_stream <- F_event$r_eq(1, mu)
    } else {
      b_stream <- F_interim$r_eq(1, lambda)
    } 
  } else {
    if (start_state) {
      b_stream <- F_event$r_gen(1, mu)
    } else {
      b_stream <- F_interim$r_gen(1, lambda)
    }    
    
  }
  
  cum_length <- b_stream 
  cum_size <- 1
  
  # add event durations and interim times until total length exceeds stream length
  while (cum_length < stream_length) {
    
    # generate random event durations and interim times
    extend_size <- ceiling(tuning * (stream_length - cum_length) / (mu + lambda))
    event_times <- F_event$r_gen(n=extend_size, mean = mu)
    interim_times <- F_interim$r_gen(n=extend_size, mean = lambda)
    
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
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
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
#'                   F_event = F_exp(), F_interim = F_exp(), 
#'                   stream_length = 100)
#'                   
#' # non-equilibrium initial conditions
#' r_behavior_stream(n = 5, mu = 3, lambda = 10,
#'                   F_event = F_gam(3), F_interim = F_gam(3),
#'                   stream_length = 100, 
#'                   equilibrium = FALSE, p0 = 0.5)

r_behavior_stream <- function(n, mu, lambda, F_event, F_interim, stream_length, 
                              equilibrium = TRUE, p0 = 0, tuning = 2) {
  
  mu_vec <- rep(mu, length.out = n)
  lambda_vec <- rep(lambda, length.out = n)
  p0_vec <- if (equilibrium) mu_vec / (mu_vec + lambda_vec) else rep(p0, length.out = n)
  
  BS <- list(stream_length = stream_length,
       b_streams = mapply(r_behavior_stream_single, 
                          mu = mu_vec,
                          lambda = lambda_vec,
                          p0 = p0_vec,
                          MoreArgs = list(F_event = F_event, F_interim = F_interim, 
                                      stream_length = stream_length, 
                                      equilibrium = equilibrium,
                                      tuning = tuning),
                          SIMPLIFY = FALSE))
  class(BS) <- "behavior_stream"
  return(BS)
}

# generate a single PIR behavior stream
r_interval_single <- function(mu, lambda, stream_length, F_event, F_interim, 
                         interval_length, coding, rest_length = 0){
  BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                          F_event = F_event, F_interim = F_interim, 
                          stream_length = stream_length)
  
  if(coding == "PIR"){
  sample <- interval_recording(BS = BS, interval_length = interval_length, rest_length = rest_length, summarize = FALSE)
  }
  
  if(coding == "WIR"){
    sample <- interval_recording(BS = BS, interval_length = interval_length, rest_length = rest_length, partial = FALSE, summarize = FALSE)
  }
  
  if(coding == "MTS"){
    sample <- momentary_time_recording(BS = BS, interval_length = interval_length, summarize = FALSE)
  }
  sample
}

#' @title Generates random partial interval recording behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as partial interval recording data with given interval length
#' and rest length.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param interval_length length of the active observation interval
#' @param rest_length length of the recording interval
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length.
#' @export
r_PIR <- function(n, mu, lambda, stream_length, F_event, F_interim, interval_length, rest_length){
  samples <- replicate(n,r_interval_single(mu = mu, lambda = lambda,
                                               stream_length = stream_length,
                                               interval_length = interval_length,
                                               rest_length = rest_length,
                                               F_event = F_event,
                                               F_interim = F_interim,
                                               coding = "PIR"))[,1,]
  t(samples)
}

#' @title Generates random whole interval recording behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as whole interval recording data with given interval length
#' and rest length.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param interval_length length of the active observation interval
#' @param rest_length length of the recording interval
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length.
#' @export

r_WIR <- function(n, mu, lambda, stream_length, F_event, F_interim, interval_length, rest_length){
  samples <- replicate(n,r_interval_single(mu = mu, lambda = lambda,
                                      stream_length = stream_length,
                                      interval_length = interval_length,
                                      rest_length = rest_length,
                                      F_event = F_event,
                                      F_interim = F_interim,
                                      coding = "WIR"))[,1,]
  t(samples)
}

#' @title Generates random momentary time sampling behavior streams
#' 
#' @description
#' Random generation of behavior streams (based on an alternating
#' renewal process) of a specified length and with specified mean event 
#' durations, mean interim times, event distribution, and interim distribution,
#' which are then coded as whole interval recording data with given interval length
#' between moments.
#' 
#' @param n number of behavior streams to generate
#' @param mu mean event duration
#' @param lambda mean interim time
#' @param F_event distribution of event durations. Must be of class \code{\link{eq_dist}}.
#' @param F_interim distribution of interim times. Must be of class \code{\link{eq_dist}}.
#' @param stream_length length of behavior stream
#' @param interval_length length of time between moments
#' @details Generates behavior streams by repeatedly drawing random event durations and 
#' random interim times from the distributions as specified, until the sum of the durations and interim
#' times exceeds the requested stream length.
#' @export
r_MTS <- function(n, mu, lambda, stream_length, F_event, F_interim, interval_length){
  samples <- replicate(n,r_interval_single(mu = mu, lambda = lambda,
                                           stream_length = stream_length,
                                           interval_length = interval_length,
                                           F_event = F_event,
                                           F_interim = F_interim,
                                           coding = "MTS"))[,1,]
  t(samples)
}