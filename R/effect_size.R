#' @title Bias corrected log ratio and variance
#' 
#' @description Estimates the bias corrected log ratio and the variance of the log ratio
#' 
#' @param observations vector of observations
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param base_level an optional value indicating the name of the baseline level. If nothing is provided, the first level in phase will be treated as the baseline level. Default is \code{NULL}
#' @param bias_correct logical value indicating if the bias-corrected log-response ratio should be used. Default is \code{null}
#' 
#' @details The \code{observations} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". 
#' If there are more than two levels in \code{phase} this function will not work. 
#' By default it is assumed that the first level in \code{phase} is the baseline level, 
#' but the baseline level can be specified via the \code{base_level} variable
#' @return A list containing two elements. The first element, \code{respRatio}, is the log-response ratio, either corrected or uncorrected depending upon the value of \code{bias_correct}. The second element, \code{variance}, is the variance of the log-response ratio.
#'
#'  @examples 
#' #Get a log ratio and variance estimate for Carl from Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#' logRespRatio(observations = outcome, phase = Phase, base_level = "No Choice"))
#' 
#' @export

logRespRatio <- function(observations, phase, base_level = NULL, bias_correct = TRUE) {
  #Data validation code follows. This is not intended to be exhaustive
  #but hopefully will catch a few of the more obvious mistakes a user might make
  if(length(which(observations < 0)) > 0 | sum(is.na(observations)) > 0) {
    stop('Values for observations must be between 0 and 1 and cannot be NA')
  }
  
  #phase_validation function is located in PIR_estimators.R
  base_level <- phase_validation(observations = observations, phase = phase, base_level = base_level)
  
  #calculates n for both samples
  base_n <- length(which(phase == base_level))
  treat_n <- length(which(phase != base_level))
  
  means <- tapply(observations, phase, fun = mean)
  variances <- tapply(observations, phase, fun = mean)
  
  if(base_level == phase[1]){
    BI = 1 #baseline index
    TI = 2 #treatment index
  } else{
    BI = 2
    TI = 1
  }

  if(means[BI] == 0){
    stop('The base mean\'s value is at the floor of 0 and this generic function does not truncate mean values')
  }
  
  if(means[TI] == 0){
    stop('The treatment mean\'s value is at the floor of 0 and this generic function does not truncate mean values')
  }

  
  if(bias_correct == TRUE){
  respRatio <- log(means[TI]) + (variances[BI]/(2 * base_n * (means[BI]^2))) 
                      - log(means[BI]) - (variances[TI]/(2 * treat_n * (means[TI]^2)))
  }else{
    respRatio <- log(means[TI]) - log(means[BI])
  }
  
  #variance of the log ratio
  variance <- (variances[BI]/(base_n * (means[BI]^2))) + 
                (variances[TI]/(treat_n * (means[TI]^2)))
  
  return(list(respRatio = respRatio, variance = variance))
}