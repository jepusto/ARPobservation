## generate behavior streams ####

#' @title An interactive alternating renewal process simulator
#' 
#' @description
#' An interactive tool that simulates single-case designs with outcomes measured by systematic direct observation.
#' The behavioral outcomes are generated from the alternating renewal process model. Both event behaviors and 
#' state behaviors are supported. Event behaviors are generated from a renewal process with gamma-distributed 
#' inter-event times. State behaviors are generated from an alternating poisson process, in which the event 
#' durations and interim times are exponentially distributed. Currently, multiple baseline and treatment reversal 
#' designs are supported. 
#' 
#' @export
#' 

ARPsimulator <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The simulator requires the shiny package. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The simulator requires the ggplot2 package. Please install it.", call. = FALSE)
  }
  
  appDir <- system.file("shiny-examples", "ARPsimulator", package = "ARPobservation")
  if (appDir == "") {
    stop("Could not find the application directory. Try re-installing ARPobservation.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}