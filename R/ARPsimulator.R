## generate behavior streams ####

#' @title An interactive alternating renewal process simulator
#'
#' @description An interactive tool that simulates single-case designs with
#' outcomes measured by systematic direct observation. The behavioral outcomes
#' are generated from the alternating renewal process model. Both event
#' behaviors and state behaviors are supported. Event behaviors are generated
#' from a renewal process with gamma-distributed inter-event times. State
#' behaviors are generated from an alternating poisson process, in which the
#' event durations and interim times are exponentially distributed. Currently,
#' multiple baseline and treatment reversal designs are supported.
#'
#' @param launch_browser Logical value indicating whether to run the tool using
#'   the system's default web browser. Defaults to \code{TRUE}.
#' @export
#' 

ARPsimulator <- function(launch_browser = TRUE) {
  
  pkgs <- c("shiny","markdown","dplyr","tidyr","ggplot2","viridis")
  loaded_pkgs <- sapply(pkgs, requireNamespace, quietly = TRUE)
  if (any(!loaded_pkgs)) {
    msg <- paste0("The simulator requires the following packages to work: ", 
                 paste(pkgs[!loaded_pkgs], collapse = ", "), 
                 ". Please install them and then try again.")
    stop(msg, call. = FALSE)
  }

  appDir <- system.file("shiny-examples", "ARPsimulator", package = "ARPobservation")
  if (appDir == "") {
    stop("Could not find the application directory. Try re-installing ARPobservation.", call. = FALSE)
  }
  
  shiny::runApp(appDir, launch.browser = launch_browser, display.mode = "normal")
}
