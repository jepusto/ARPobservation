#--------------------------------------------
# Plot method for behavior_stream objects
#--------------------------------------------

get_segments <- function(bs, stream_length) {
  transitions <- c(0, bs$b_stream, stream_length)
  s <- trunc((length(transitions) - 1) / 2)
  start <- transitions[2 * (1:s) - bs$start_state]
  end <- transitions[2 * (1:s) - bs$start_state + 1]
  data.frame(stream = bs$index, start, end)
}

#' @title Plot method for \code{behavior_stream} objects
#'   
#' @description Creates a graphical representation of a set of simulated
#' behavior streams.
#' 
#' @param x object of class \code{behavior_stream}
#'   
#' @details The plot is created using \code{ggplot} from the ggplot2 package, which must be installed.
#' @export
#' 
#' @return An object of class \code{ggplot}.
#'   
#' @examples
#' b_streams <- r_behavior_stream(n = 5, mu = 3, lambda = 10, 
#'                                F_event = F_exp(), F_interim = F_exp(), 
#'                                stream_length = 100)
#' plot(b_streams)

plot.behavior_stream <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This function needs the ggplot2 package to work. Please install it.", call. = FALSE)
  } else {
    require(ggplot2, quietly = TRUE)
  }
 
  streams <- length(x$b_streams)
  stream_dat <- data.frame(stream = 1:streams, start = 0, end = x$stream_length) 
  
  for (i in 1:streams) x$b_streams[[i]]$index <- i
  segment_dat <- lapply(x$b_streams, get_segments, stream_length = x$stream_length)
  for (i in 1:streams) {
    print(x$b_streams[[i]])
    print(segment_dat[[i]])
  }
  segment_dat <- do.call(rbind, segment_dat)
  
  ggplot(stream_dat, aes(x = start, xend = end, y = stream, yend = stream)) +
          geom_segment() + 
          geom_segment(data = segment_dat, color = "blue", size = 2) + 
          labs(x = "session time", y = "stream") + 
          theme_minimal()
}