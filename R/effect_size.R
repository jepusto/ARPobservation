#' @title Bias corrected log ratio and variance
#' 
#' @description Estimates the bias corrected log ratio and the variance of the log ratio
#' 
#' @param observations vector of observations
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param first_level_base Logical value indicating if the first level of \code{phase} is the baseline level. Default is \code{TRUE}.
#' 
#' @details The \code{observations} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". 
#' If there are more than two levels in \code{phase} this function will not work. 
#' By default it is assumed that the first level in \code{phase} is the baseline level, 
#' but this can be changed by setting \code{first_level_base} to \code{FALSE}.
#' @return A list containing two elements. The first element, \code{bias_corrected_R}, is the bias-corrected log ratio. The second element, \code{variance_R}, is the variance of the bias-corrected log ratio.
#'
#'  @examples 
#' #Get a log ratio and variance estimate for Carl from Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#'  bc_log_ratio(observations = outcome, phase = Phase,
#'  first_level_base = FALSE))
#' 
#' @export

bc_log_ratio <- function(observations, phase, first_level_base = TRUE) {
  #Data validation code follows. This is not intended to be exhaustive
  #but hopefully will catch a few of the more obvious mistakes a user might make
  if(length(which(observations < 0)) > 0 | sum(is.na(observations)) > 0) {
    stop('Values for observations must be between 0 and 1 and cannot be NA')
  }
  
  #recasts phase as a factor if it is just a vector
  phase <- factor(phase)
  
  #Checks to be sure there are only two levels
  if(length(levels(phase)) != 2){
    stop('This function requires exactly two levels in the phase variable')
  }
  
  #Sets a value for the name of the level of the base phase
  if(first_level_base == TRUE){
    base_phase <- levels(phase)[which(levels(phase) == phase[1])]
  } else{
    base_phase <- levels(phase)[which(levels(phase) != phase[1])]
  }
  
  #calculates n for both samples
  base_n <- length(which(phase == base_phase))
  treat_n <- length(which(phase != base_phase))
  
  #calculate mean and variance for the first sample
  base_mean <- mean(observations[which(phase == base_phase)])
  base_variance <- var(observations[which(phase == base_phase)])
  
  #This checks to see if the mean and variance reach the floor
  if(base_mean == 0){
    stop('The base mean\'s value is at the floor of 0 and this generic function does not truncate mean values')
  }
  
  treat_mean <- mean(observations[which(phase != base_phase)])
  treat_variance <- var(observations[which(phase != base_phase)])
  
  if(treat_mean == 0){
    stop('The treatment mean\'s value is at the floor of 0 and this generic function does not truncate mean values')
  }
  
  #Natural log of the ratio of the two means
  bias_corrected_R <- log(treat_mean) + (base_variance/(2 * base_n * (base_mean^2))) 
                      - log(base_mean) - (treat_variance/(2 * treat_n * (treat_mean^2)))
  
  
  #variance of the log ratio
  variance_R <- (base_variance/(base_n * (base_mean^2))) + (treat_variance/(treat_n * (treat_mean^2)))
  
  return(list(bias_corrected_R = bias_corrected_R, variance_R = variance_R))
}