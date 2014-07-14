#Many of the functions in this code re-use code from the previous functions. Commenting
#is not duplicated across functions so if something appears to be insufficientily commented
#check previous functions for explanations.


#' @title Prevalence bounds and confidence interval
#' 
#' @description Calculates a bound for the log of the prevalence ratio of two samples (referred to as baseline and treatment)
#' based on partial interval recording (PIR) data.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param mu_L the lower limit on the mean event duration
#' @param active_length length of the active observation interval
#' @param intervals the number of intervals in the sample of observations. Default is \code{NA}.
#' @param first_level_base Logical value indicating if the first level of \code{phase} is the baseline level. Default is \code{TRUE}.
#' @param conf_level Desired coverage rate of the calculated confidence interval. Default is \code{.95}.
#' @param exponentiate Logical value indicating if the log of the bounds and the confidence interval should be exponentiated. Default is \code{FALSE}.
#' 
#' @details The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". 
#' If there are more than two levels in \code{phase} this function will not work. 
#' By default it is assumed that the first level in \code{phase} is the baseline level, 
#' but this can be changed by \code{first_level_base} to \code{FALSE}.
#' 
#' \code{mu_L} is the lower limit on the mean event durations. This is a single value assumed to hold for both samples of behavior
#' 
#' \code{active_length} This is the total active observation length. If the intervals are 15 seconds long but 5 seconds of each interval is reserved for recording purposes, \code{active_length= 10}. Times are often in seconds, but can be in any time unit.
#' 
#' \code{intervals} is the number of intervals in the observations. This is a single value and is assumed to be constant across both samples and all observations. This value is only relevant if the mean of one of the samples is at the floor or ceiling of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list containing four elements
#'
#'  @examples 
#' #Get an estimate and CI for Carl from Moes dataset
#' data(FMexample)
#' with(subset(Moes, Case == "Carl"),
#'  prevalence_bounds(PIR = outcome, phase = Phase, mu_L = 10, active_length = 10, first_level_base = FALSE))
#' 
#' @export

prevalence_bounds <- function(PIR, phase, mu_L, active_length, intervals = NA,
                              first_level_base = TRUE,
                              conf_level = 0.95, exponentiate = FALSE) {
  #Data validation code follows. This is not intended to be exhaustive
  #but hopefully will catch a few of the more obvious mistakes a user might make
  
  #Checks to be sure that the values for PIR are between 0 and 1, inclusive
  #and none of the values are NA
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
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
  base_mean <- mean(PIR[which(phase == base_phase)])
  base_variance <- var(PIR[which(phase == base_phase)])
  
  #This checks to see if the mean and variance reach the floor or ceiling.
  #If not value for intervals is provided (making truncation impossible) the function
  #stops. Otherwise it will truncate the values appropriately
  if((base_mean == 0 | base_mean == 1) & is.na(intervals) == TRUE){
    stop('The baseline mean PIR value is at the floor of 0 or ceiling of 1 and no value of K has been provided to truncate the means.')
    }
  if(base_mean == 0){
    base_mean <- 1/(base_n * intervals)
  }   
  if(base_mean == 1){
    base_mean <- 1-(1/(base_n * intervals))
  }
  
  treat_mean <- mean(PIR[which(phase != base_phase)])
  treat_variance <- var(PIR[which(phase != base_phase)])
  
  if((treat_mean == 0 | treat_mean == 1) & is.na(intervals) == TRUE){
    stop('The treatment mean PIR value is at the floor of 0 or ceiling of 1 and no value for K has been provided to truncate the means.')
  }
  if(treat_mean == 0){
    treat_mean <- 1/(treat_n * intervals)
  }   
  if(treat_mean == 1){
    treat_mean <- 1-(1/(treat_n * intervals))
  }
  
  #Natural log of the ratio of the two means
  R <- log(treat_mean) - log(base_mean)
  
  h <- log(mu_L + active_length) - log(mu_L)
  
  #calculate lower and upper bounds
  lower_bound <- R - h
  upper_bound <- R + h
  
  #variance of the log ratio
  variance_R <- (base_variance/(base_n * (base_mean^2))) + (treat_variance/(treat_n * (treat_mean^2)))
  
  #calculate the CI
  lower_CI <- R - (h + (qnorm(1-(1-conf_level)/2)*sqrt(variance_R)))
  upper_CI <- R + (h + (qnorm(1-(1-conf_level)/2)*sqrt(variance_R)))
  
  #exponentiates the values, if desired
  if(exponentiate == TRUE){
    lower_bound <- exp(lower_bound)
    upper_bound <- exp(upper_bound)
    lower_CI <- exp(lower_CI)
    upper_CI <- exp(upper_CI)
  }
  
  return(list(lower_bound = lower_bound, upper_bound = upper_bound, 
              lower_CI = lower_CI, upper_CI = upper_CI))
}

#' Incidence bounds and confidence interval
#' 
#' @description This function calculates a bound for the log of the incidence ratio of two samples (referred to as baseline and treatment) of partial interval recording (PIR) data as well as an approximate confidence interval
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param mu_U the upper limit on the mean event duration
#' @param p the probability that the interim time between behavioral events is less than the active interval
#' @param active_length length of the active observation interval
#' @param intervals \code{NA} (default) or the number of intervals in the sample of observations
#' @param first_level_base logical. Is the first level in the factor or vector \code{phase} the baseline level?
#' @param conf_level \code{.95} (default) or the desired nominal confidence level of the calculated confidence interval
#' @param exponentiate logical. Should the log of the bounds and the confidence interval be exponentiated?
#' 
#' @details The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". If there are more than two levels in \code{phase} this function will not work. By default it is assumed that the first level in \code{phase} is the baseline level. This behavior can be changed by \code{first_level_base} to \code{FALSE}.
#' 
#' \code{mu_U} is the upper limit on the mean event durations. This is a single value assumed to hold for both samples of behavior
#' 
#' \code{active_length} This is the total active observation length. If the intervals are 15 seconds long but 5 seconds of each interval is reserved for recording purposes, \code{active_length= 10}. Times are often in seconds, but can be in any time unit.
#' 
#' \code{intervals} is the number of intervals in the observations. This is a single value and is assumed to be constant across both samples and all observations. This value is only relevant if the mean of one of the samples is at the floor or ceiling of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list containing four elements
#' 
#' @examples 
#' #get an estimate and CI for Ahmad from the Dunlap dataset
#' data(FMexample)
#' PIR <- Dunlap[1:16, 4]
#' phase <- Dunlap[1:16, 2]
#' incidence_bounds(PIR = PIR, phase = phase, mu_U = 10, p = .15, active_length = 10, intervals = 60)
#' 
#' @export
incidence_bounds <- function(PIR, phase, mu_U, p, active_length, intervals = NA, 
                             first_level_base = TRUE,
                             conf_level = 0.95, exponentiate = FALSE) {
  
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  phase <- factor(phase)
  
  if(length(levels(phase)) != 2){
    stop('This function requires exactly two levels in the phase variable')
  }
  
  if(first_level_base == TRUE){
    base_phase <- levels(phase)[which(levels(phase) == phase[1])]
  } else{
    base_phase <- levels(phase)[which(levels(phase) != phase[1])]
  }
  
  base_n <- length(which(phase == base_phase))
  treat_n <- length(which(phase != base_phase))
  
  base_mean <- mean(PIR[which(phase == base_phase)])
  base_variance <- var(PIR[which(phase == base_phase)])
  
  if((base_mean == 0 | base_mean == 1) & is.na(intervals) == TRUE){
    stop('The baseline mean PIR value is at the floor of 0 or ceiling of 1 and no value of K has been provided to truncate the means.')
  }
  if(base_mean == 0){
    base_mean <- 1/(base_n * intervals)
  }   
  if(base_mean == 1){
    base_mean <- 1-(1/(base_n * intervals))
  }
  
  treat_mean <- mean(PIR[which(phase != base_phase)])
  treat_variance <- var(PIR[which(phase != base_phase)])
  
  if((treat_mean == 0 | treat_mean == 1) & is.na(intervals) == TRUE){
    stop('The treatment mean PIR value is at the floor of 0 or ceiling of 1 and no value for K has been provided to truncate the means.')
  }
  if(treat_mean == 0){
    treat_mean <- 1/(treat_n * intervals)
  }   
  if(treat_mean == 1){
    treat_mean <- 1-(1/(treat_n * intervals))
  }
  
  R <- log(treat_mean) - log(base_mean)
  
  h <- log(mu_U + active_length) - log(1-p) - log(active_length)
  
  lower_bound <- R - h
  upper_bound <- R + h
  
  variance_R <- (base_variance/(base_n * (base_mean^2))) + (treat_variance/(treat_n * (treat_mean^2)))
  
  lower_CI <- R - (h + (qnorm(1-(1-conf_level)/2) * sqrt(variance_R)))
  upper_CI <- R + (h + (qnorm(1-(1-conf_level)/2) * sqrt(variance_R)))
 
  if(exponentiate == TRUE){
    exp(lower_bound)
    exp(upper_bound)
    exp(lower_CI)
    exp(upper_CI)
  }
                   
  return(list(lower_bound = lower_bound, upper_bound = upper_bound, 
              lower_CI = lower_CI, upper_CI = upper_CI))
}

#' @title Interim bounds and confidence interval
#' 
#' @description This function calculates a bound for the log of the ratio of interim time of two samples (referred to as baseline and treatment) of partial interval recording (PIR) data as well as an approximate confidence interval
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param first_level_base logical. Is the first level in the factor or vector \code{phase} the baseline level?
#' @param conf_level \code{.95} (default) or the desired nominal confidence level of the calculated confidence interval
#' @param intervals \code{NA} (default) or the number of intervals in the sample of observations
#' @param exponentiate logical. Should the log of the bounds and the confidence interval be exponentiated?
#' @details The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". If there are more than two levels in \code{phase} this function will not work. By default it is assumed that the first level in \code{phase} is the baseline level. This behavior can be changed by \code{first_level_base} to \code{FALSE}.
#' \code{intervals} is the number of intervals in the observations. This is a single value and is assumed to be constant across both samples and all observations. This value is only relevant if the mean of one of the samples is at the floor or ceiling of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list containing four elements
#' 
#' @examples 
#' #get an estimate and CI for Carl from the Moes dataset
#' data(FMexample)
#' PIR <- Moes[1:20,4]
#' phase <- Moes[1:20, 2]
#' interim_bounds(PIR = PIR, phase = phase, first_level_base = FALSE)
#' 
#' @export
interim_bounds <- function(PIR, phase, first_level_base = TRUE, 
                           conf_level = 0.95, intervals = NA, exponentiate = FALSE) {
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  phase <- factor(phase)
  
  if(length(levels(phase)) != 2){
    stop('This function requires exactly two levels in the phase variable')
  }
  
  if(first_level_base == TRUE){
    base_phase <- levels(phase)[which(levels(phase) == phase[1])]
  } else{
    base_phase <- levels(phase)[which(levels(phase) != phase[1])]
  }
  
  base_mean <- mean(PIR[which(phase == base_phase)])
  base_variance <- var(PIR[which(phase == base_phase)])
  
  if((base_mean == 0 | base_mean == 1) & is.na(intervals) == TRUE){
    stop('The baseline mean PIR value is at the floor of 0 or ceiling of 1 and no value of K has been provided to truncate the means.')
  }
  if(base_mean == 0){
    base_mean <- 1/(base_n * intervals)
  }   
  if(base_mean == 1){
    base_mean <- 1-(1/(base_n * intervals))
  }
  
  treat_mean <- mean(PIR[which(phase != base_phase)])
  treat_variance <- var(PIR[which(phase != base_phase)])
  
  if((treat_mean == 0 | treat_mean == 1) & is.na(intervals) == TRUE){
    stop('The treatment mean PIR value is at the floor of 0 or ceiling of 1 and no value for K has been provided to truncate the means.')
  }
  if(treat_mean == 0){
    treat_mean <- 1/(treat_n * intervals)
  }   
  if(treat_mean == 1){
    treat_mean <- 1-(1/(treat_n * intervals))
  }
  
  base_n <- length(which(phase == base_phase))
  treat_n <- length(which(phase != base_phase))
  
  #the logit and complimentary log-log functions
  logit <- function(x) {log(x) - log(1-x)}
  cll <- function(x) {log(-1 * log(1-x))}
  
  #This code checks to see which transformation is the appropriate value for
  #the lower and upper bounds
  if(base_mean > treat_mean){
    f_lower <- logit(treat_mean) - logit(base_mean)
    f_upper <- cll(treat_mean) - cll(base_mean)
  }else{
    f_lower <- cll(treat_mean) - cll(base_mean)
    f_upper <- logit(treat_mean) - logit(base_mean)
  }
  
  
  #This is the variance of the log ratio of the sample means
  var_LOR <- (base_variance/(base_n * (base_mean^2) * (1- base_mean)^2)) + 
             (treat_variance/(treat_n * (treat_mean^2) * (1- treat_mean)^2))
  
  #This is the variance of the complimentary log-log ratio
  var_CLR <- (base_variance/(base_n * (1- base_mean)^2 * (log(1-base_mean))^2)) + 
             (treat_variance/(treat_n * (1- treat_mean)^2 *(log(1-treat_mean))^2))
  
  #The ratio and z_conf is used in determining which variance is appropriate
  cll_ratio <- cll(treat_mean) - cll(base_mean)
  
  z_conf <- qnorm(1-(1-conf_level)/2)
  
  if(cll_ratio <= (z_conf * sqrt(var_LOR))){
    var_f_lower <- var_LOR
  } else{
    var_f_lower <- var_CLR
  }
  
  if(cll_ratio < -(z_conf * sqrt(var_LOR))){
    var_f_upper <-  var_CLR
  }else{
    var_f_upper <- var_LOR
  }
  
  lower_CI <- f_lower - (z_conf * sqrt(var_f_lower))
  upper_CI <- f_upper + (z_conf * sqrt(var_f_upper))
  
  if(exponentiate == TRUE)
  {
    f_lower <- exp(f_lower)
    f_upper <- exp(f_upper)
    lower_CI <- exp(lower_CI)
    upper_CI <- exp(upper_CI)
  }
  
  return(list(lower_bound = f_lower, upper_bound = f_upper,
              lower_CI = lower_CI, upper_CI = upper_CI))
}

#estimate the value of zeta based on the expected value of zeta and the estimate of phi
#and the sample mean
PIR_Zeta <- function(phi, ExY, active) -1 * (1-phi) * log((ExY-1)/(phi-1)) / active

#estimate the variance of the sample mean using the expected value of the variance
#and the estimates of phi and zeta
PIR_VarEx <- function(phi, ExY, active, L, K){
  
  zeta = PIR_Zeta(phi, ExY, active)
  
  if (phi == 0) {
    vary = ExY * (1- ExY) / K 
  } else {
    k = 1:(K-1) #vector for summation term
    vary = ((1/K) * ExY * (1-ExY)) * ( 1 + (2 * phi / K / ExY) * sum((K - k) * exp((zeta * active / phi)  - zeta * k * L /phi / (1- phi) / K)) )
  }  
  
  return(vary)
}

#Finds the root for phi using the variance and mean of the sample
PIR_Phi <- function(ExY, VarY, nObs, active, L, K) {
  fun <- function(phi) PIR_VarEx(phi, ExY, active, L, K) / VarY - 1
  uniroot(fun, interval = c(0, ExY), tol = .Machine$double.eps^0.5)$root
}

#generates multiple samples of PIR data using the r_behavior_stream and interval_recording
#functions from ARPobservation. Useful for simulating or performing bootstraps
generatePIRData <- function(nObs, phi, zeta, active, K, rest = 0, iterations = 1) {
  
  #necessary to recast as numbers because PIR_MOM sends them as lists with one element
  phi <- as.numeric(phi) 
  zeta <- as.numeric(zeta)
  #calculating the means of the distributions of event durations and interim
  #times from estimates of prevalence and incidence
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  
  #the full interval length is the sum of the active and rest (recording) lengths
  intervalLength = active + rest
  
  #this block generates a matrix of PIR data with each column corresponding to a
  #a sample of n = nObs observations
  sampleObs <- replicate(iterations, {
    BS <- r_behavior_stream(n = nObs, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = intervalLength * K)
    interval_recording(BS = BS, interval_length = intervalLength, rest_length = rest)
  })
  
  mean <- colMeans(sampleObs)
  
  variance <- apply(sampleObs, 2, var)
  
  sampleData <- data.frame(mean = mean, variance = variance)
  
  return(sampleData)
}

#A packaging function that invokes the previous set of functions to neatly
#provide an estimate of phi and zeta for a sample or vector of samples
PIR_inv <- function(ExY, VarY, nObs, active, L, K, varTrunc = 1 / (nObs * K^2)) {
  
  ExY <- pmin(pmax(ExY, 1 / (nObs * K)), 1 - 1 / (nObs * K))
  VarY <- pmin(pmax(VarY, ExY * (1 - ExY) / K + varTrunc), ExY * (1 - ExY) - varTrunc)
  
  phi <- mapply(PIR_Phi, ExY = ExY, VarY = VarY, 
                nObs = nObs, active = active, L = L, K = K)
  
  zeta <- mapply(PIR_Zeta, phi = phi, ExY = ExY, active = active)
  
  ests <- data.frame(phi = phi, zeta = zeta)
  
  return (ests)
  
}

#bootstraps the confidence intervals for a pair of samples with estimates of phi and zeta
PIRbootstrappair<- function(nObs, phi, zeta, active, rest, K, iterations, alpha,
                            seed = NULL){
  #set seed if one is supplied
  if (!is.null(seed)) set.seed(seed)

  #the total length of the observation is the total interval length times the number of intervals
  L = (active + rest) * K
  
  #Simulate first set of parms
  sampleData0 <- generatePIRData(nObs = nObs[1], phi = phi[1], zeta = zeta[1],
                                 active = active, K = K, rest = rest, 
                                 iterations = iterations)
  
  #This trims the means so that they are not at 0 or 1
  sampleData0$meanTrim <- pmin(pmax(sampleData0$mean, 1 / (nObs[1] * K)), 
                               1 - 1 / (nObs[1] * K))
  
  #Returns a dataframe of estimates of phi and zeta
  ests0 <- with(sampleData0, PIR_inv(ExY = meanTrim, VarY = variance, 
                                     nObs = nObs[1], active = active, K = K, 
                                     L = L))
  #truncates the estimates so they are not at zero
  ests0$ptrunc <- pmax(ests0$phi, 1/(nObs*K))
  ests0$ztrunc <- pmax(ests0$zeta, 1/(active*nObs*K))
  
  #get the confidence interval bounds for the first sample
  pbounds0 <- quantile(ests0$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds0 <- quantile(ests0$zeta, probs = c(alpha/2, 1-alpha/2))
  
  #simulate second set of parms
  sampleData1 <- generatePIRData(nObs = nObs[2], phi = phi[2], zeta = zeta[2],
                                 active = active, K = K, rest = rest, 
                                 iterations = iterations)
  
  sampleData1$meanTrim <- pmin(pmax(sampleData1$mean, 1 / (nObs[2] * K)), 
                               1 - 1 / (nObs[2] * K))
  
  ests1 <- with(sampleData1, PIR_inv(ExY = meanTrim, VarY = variance, 
                                     nObs = nObs[2], active = active, K = K, 
                                     L = L))
  ests1$ptrunc <- pmax(ests1$phi, 1/(nObs*K))
  ests1$ztrunc <- pmax(ests1$zeta, 1/(active*nObs*K))
  pbounds1 <- quantile(ests1$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds1 <- quantile(ests1$zeta, probs = c(alpha/2, 1-alpha/2))
  
  #calculate log ratio value and confidence interval bounds
  plogratio <- log(ests1$ptrunc/ests0$ptrunc)
  plogbounds <- quantile(plogratio, probs = c(alpha/2, 1-alpha/2))
  
  zlogratio <- log(ests1$ztrunc/ests0$ztrunc)
  zlogbounds <- quantile(zlogratio, probs = c(alpha/2, 1-alpha/2))
  
  results <- cbind(phi = c(phi[1], phi[2], log(as.numeric(phi[2])/as.numeric(phi[1]))),
                   phi_lower_CI = c(pbounds0[1], pbounds1[1], plogbounds[1]),
                   phi_upper_CI = c(pbounds0[2], pbounds1[2], plogbounds[2]),
                   zeta = c(zeta[1], zeta[2], log(as.numeric(zeta[2])/as.numeric(zeta[1]))),
                   zeta_lower_CI = c(zbounds0[1], zbounds1[1], zlogbounds[1]),
                   zeta_upper_CI = c(zbounds0[2], zbounds1[2], zlogbounds[2]))
  
  return(results)
  
}

#' @title Moment estimator for prevalence and incidence and bootstrap confidence intervals
#' 
#' @description This function provides estimates of prevalance and incidence for two samples as well as the log ratios. It provides boostrap confidence intervals using the functions \code{\link{r_behavior_stream}} and \code{\link{interval_recording}} to generate the replicates.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param intervals the number of intervals in the sample of observations
#' @param interval_length the total length of the intervals
#' @param rest_length \code{0} (default) or length of the portion of the interval devoted to recording
#' @param first_level_base logical. Is the first level in the factor or vector \code{phase} the baseline level?
#' @param Bootstraps \code{2000} (default) or the desired number of bootstrap replicates
#' @param conf_level \code{.95} (default) or the desired nominal confidence level of the calculated confidence interval
#' @param seed \code{NULL} (default) or a seed value set in order to make bootstrap results reproducible
#' 
#' @details The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". If there are more than two levels in \code{phase} this function will not work. By default it is assumed that the first level in \code{phase} is the baseline level. This behavior can be changed by \code{first_level_base} to \code{FALSE}.
#' 
#' \code{intervals}, \code{interval_length}, and \code{rest_length} are all single values that are assumed to be held constant across both samples and all observation sessions.
#' 
#' \code{interval_length} This is the length of each individual interval. Sometimes a portion of the interval is set aside for recording purposes. The length of that recording interval is what the value of \code{rest_length} should be set to, although the default assumption is that there is no recording time. The length of the active interval is calculated to be \code{interval_length - rest_length}.
#' 
#' \code{bootstraps} defaults to 2000. This is believed to be a sufficient number. At the default, PIR_MOM takes just under six seconds to run on an Intel Core i5-2410M processor.
#' 
#' \code{seed} defaults to \code{NULL}, but allows for reproducible code.
#' 
#' @return A dataframe with three rows corresponding to baseline, treatment, and treatment/baseline and six columns
#' 
#' @examples 
#' #get estimate and CIs for Carl from the Moes dataset
#' data(FMexample)
#' PIR <- Moes[1:20,4]
#' phase <- Moes[1:20, 2]
# 'PIR_MOM(PIR = PIR, phase = phase, intervals = 80, interval_length = 15, rest_length = 5,first_level_base= FALSE, seed = 149568373)
#' 
#' @export
PIR_MOM <- function(PIR, phase, intervals, interval_length, rest_length = 0, first_level_base = TRUE, Bootstraps = 2000, conf_level = 0.95, seed = NULL) {
  
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  phase <- factor(phase)
  
  if(length(levels(phase)) != 2){
    stop('This function requires exactly two levels in the phase variable')
  }
  
  if(first_level_base == TRUE){
    base_phase <- levels(phase)[which(levels(phase) == phase[1])]
  } else{
    base_phase <- levels(phase)[which(levels(phase) != phase[1])]
  }
  
  base_mean <- mean(PIR[which(phase == base_phase)])
  base_variance <- var(PIR[which(phase == base_phase)])
  
  if((base_mean == 0 | base_mean == 1) & is.na(intervals) == TRUE){
    stop('The baseline mean PIR value is at the floor of 0 or ceiling of 1 and no value of K has been provided to truncate the means.')
  }
  if(base_mean == 0){
    base_mean <- 1/(base_n * intervals)
  }   
  if(base_mean == 1){
    base_mean <- 1-(1/(base_n * intervals))
  }
  
  treat_mean <- mean(PIR[which(phase != base_phase)])
  treat_variance <- var(PIR[which(phase != base_phase)])
  
  if((treat_mean == 0 | treat_mean == 1) & is.na(intervals) == TRUE){
    stop('The treatment mean PIR value is at the floor of 0 or ceiling of 1 and no value for K has been provided to truncate the means.')
  }
  if(treat_mean == 0){
    treat_mean <- 1/(treat_n * intervals)
  }   
  if(treat_mean == 1){
    treat_mean <- 1-(1/(treat_n * intervals))
  }
  
  base_n <- length(which(phase == base_phase))
  treat_n <- length(which(phase != base_phase))
  
  #Get an estimate of phi and zeta using the moment estimator for baseline and treatment phases
  base_ests <- PIR_inv(ExY = base_mean, VarY = base_variance,
                       nObs = base_n, active = interval_length - rest_length,
                       L = intervals * interval_length, K = intervals)  
  

  treat_ests <- PIR_inv(ExY = treat_mean, VarY = treat_variance,
                      nObs = treat_n, active = interval_length - rest_length, 
                      L = intervals * interval_length, K = intervals)
  
  #Bootstrap the confidence intervals
  results <- PIRbootstrappair(nObs = c(base_n, treat_n),
                   phi = c(base_ests[1], treat_ests[1]),
                   zeta = c(base_ests[2], treat_ests[2]),
                   active = interval_length - rest_length,
                   rest = rest_length,
                   K = intervals,
                   iterations = Bootstraps,
                   alpha = 1-conf_level,
                   seed = seed)
  
  #Set the row names appropriately based on the levels in the "phase" variable
  row.names(results) <- c(base_phase, 
                          levels(phase)[which(levels(phase) != base_phase)], 
                          paste(levels(phase)[which(levels(phase) != base_phase)], base_phase, sep = "/"))
  
  return(results)
}

#' Dunlap et al.(1994) and Moes(1998) data
#' 
#' Single case design data from two different studies examining the impact of choice-making on disruptive behavior in academic settings. Thus, baseline is "no choice" and treatment was "choice." Data were extracted from the figures in the publicatons.
#' 
#' For the Dunlap dataset Ahmad and Sven were measured with an active interval of 10s with 5s for recording while Wendall had an active interval length of 15s with no time for recording. Sessions were 15 minutes long and each summary measurement was based on K = 60 intervals.
#' 
#' For the Moes dataset all participants were observed with an active interval of 10s with 5s. Sessions were 20 minutes long and each summary measurement was based on K = 80 intervals.
#' 
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 53940 rows and 10 variables
#' @name FMexample
#' @aliases Moes Dunlap
#' @references 
#' Dunlap, G., DePerczel, M., Clarke, S., Wilson, D., Wright, S., White, R., & Gomez, A. (1994). Choice making to promote adaptive behavior for students with emotional and behavioral challenges. Journal of Applied Behavior Analysis, 27 (3), 505-518.
#' 
#' Moes, D. R. (1998). Integrating choice-making opportunities within teacher-assigned academic tasks to facilitate the performance of children with autism. Research and Practice for Persons with Severe Disabilities, 23 (4), 319-328.
NULL
