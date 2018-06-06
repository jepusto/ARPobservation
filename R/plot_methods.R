#--------------------------------------------
# Plot method for behavior_stream objects
#--------------------------------------------

get_segments <- function(bs, stream_length) {
  transitions <- c(0, bs$b_stream, stream_length)
  s <- trunc((length(bs$b_stream) + bs$start_state + 1) / 2)
  start_time <- transitions[2 * (1:s) - bs$start_state]
  end_time <- transitions[2 * (1:s) - bs$start_state + 1]
  data.frame(stream = bs$index, start_time, end_time)
}

#' @title Plot method for \code{behavior_stream} objects
#'   
#' @description Creates a graphical representation of a set of simulated 
#'   behavior streams.
#'   
#' @param x object of class \code{behavior_stream}
#' @param session_color character string indicating the color of the lines that
#'   represent session time. Default is black.
#' @param episode_color character string indicating the color of the bars that
#'   represent episode durations. Default is blue.
#' @param episode_thickness numeric value indicating the thickness of the bars
#'   that represent epsiode durations. Default is 2.
#' @param ... Further arguments, not used for this method.
#'   
#' @details The plot is created using \code{ggplot} from the ggplot2 package, 
#'   which must be installed.
#' @export
#' 
#' @return An object of class \code{ggplot}.
#'   
#' @examples
#' b_streams <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                                F_event = F_exp(), F_interim = F_exp(), 
#'                                stream_length = 100)
#' plot(b_streams)

plot.behavior_stream <- function(x, session_color = "black", episode_color = "blue", episode_thickness = 2, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This function needs the ggplot2 package to work. Please install it.", call. = FALSE)
  } 

  streams <- length(x$b_streams)
  stream_dat <- data.frame(stream = 1:streams, start_time = 0, end_time = x$stream_length) 
  
  for (i in 1:streams) x$b_streams[[i]]$index <- i
  segment_dat <- lapply(x$b_streams, get_segments, stream_length = x$stream_length)
  segment_dat <- do.call(rbind, segment_dat)
  
  ggplot2::ggplot(stream_dat, ggplot2::aes(x = start_time, xend = end_time, y = stream, yend = stream)) +
    ggplot2::geom_segment(color = session_color) + 
    ggplot2::scale_y_discrete(breaks = 1:streams) + 
    ggplot2::coord_cartesian(xlim = c(0, x$stream_length)) + 
    ggplot2::geom_segment(data = segment_dat, color = episode_color, size = episode_thickness) + 
    ggplot2::labs(x = "session time", y = "stream") + 
    ggplot2::theme_minimal()
}
