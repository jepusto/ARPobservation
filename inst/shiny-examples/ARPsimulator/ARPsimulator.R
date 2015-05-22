#------------------------------
# Design matrix
#------------------------------

phase_design <- function(design, cases, phase_pairs, sessions_TR, 
                         sessions_MB, phase_change_list, samples) {
  
  case_names <- paste("Case", LETTERS[1:cases])
  
  if (design=="Treatment Reversal") {
    sample_nr <- 1:samples
    session_nr <- 1:(phase_pairs * 2 * sessions_TR)
    phase_changes <- sessions_TR * (1:(2 * phase_pairs - 1))
    phase <- rep(1:(2 * phase_pairs), each = sessions_TR)
    trt <- rep(rep(c("Base","Treat"), each = sessions_TR), times = phase_pairs)
  } else {
    sample_nr <- 1:samples
    session_nr <- 1:sessions_MB
    phase_changes <- as.numeric(strsplit(phase_change_list, split = ",")[[1]])
    phase_changes <- rep(phase_changes, length.out = cases)
    phase <- unlist(lapply(phase_changes, function(t) rep(1:2, times = c(t, sessions_MB - t))))
    trt <- c("Base","Treat")[phase]
  }
  
  dat <- expand.grid(session = session_nr, case = case_names, sample = sample_nr)
  dat$phase <- factor(phase)
  dat$trt <- factor(trt)
  
  return(list(dat = dat, phase_changes = phase_changes))
}

#----------------------------
# Generate outcome data
#----------------------------

go <- function(behavior, freq, duration, interim_time) {
  (behavior == "Event behavior" & !is.na(freq)) | 
       (behavior == "State behavior" & !is.na(duration) & !is.na(interim_time))
}

simulate_measurements <- function(trt, behavior, freq, dispersion, freq_change, 
                                  duration, interim_time, duration_change, interim_change,
                                  system, interval_length, session_length) {
  
  if (behavior == "Event behavior") {
    mu <- c(0,0)
    lambda <- 1 / (freq * c(1, freq_change / 100 + 1))
    F_event <- F_const()
    F_interim <- F_gam(shape = 1 / dispersion)
  } else {
    mu <- duration * c(1, duration_change / 100 + 1) / 60
    lambda <- interim_time * c(1, interim_change / 100 + 1) / 60
    F_event <- F_exp()
    F_interim <- F_exp()
  }
  
  BS <- r_behavior_stream(n = length(trt), mu = mu[trt], lambda = lambda[trt],
                          F_event = F_event, F_interim = F_interim, stream_length = session_length)
  
  Y <- switch(system,
              "Frequency counting" = event_counting(BS),
              "Continuous recording" = 100 * continuous_duration_recording(BS),
              "Momentary time sampling" = 100 * momentary_time_recording(BS, interval_length / 60),
              "Partial interval recording" = 100 * interval_recording(BS, interval_length / 60, partial = TRUE),
              "Whole interval recording" = 100 * interval_recording(BS, interval_length / 60, partial = FALSE)
  )
  return(Y)
}


#----------------------------
# Create SCD graph
#----------------------------

graph_SCD <- function(dat, design, phase_changes, system) {
  
  Y_lab <- switch(system,
                  "Frequency counting" = "Events",
                  "Continuous recording" = "% Session time",
                  "Momentary time sampling" = "% Moments",
                  "Partial interval recording" = "% Partial intervals",
                  "Whole interval recording" = "% Whole intervals"
  )
  Y_range <- c(0, ifelse(system=="Frequency counting", max(dat$Y), 100))
  
  SCD_graph <- ggplot(dat, aes(session, Y, color = trt, group = interaction(phase, sample))) + 
    geom_point() + 
    geom_line() + 
    facet_grid(case ~ .) + 
    coord_cartesian(ylim = Y_range) + 
    theme_bw() + theme(legend.position = "bottom") + 
    labs(color = "Phase", y = Y_lab) 
  
  if (design == "Multiple Baseline") {
    phase_dat <- data.frame(case = levels(dat$case), phase_change = phase_changes + 0.5)
    phase_change_lines <- geom_vline(data = phase_dat, aes(xintercept = phase_change), linetype = "dashed")
  } else {
    phase_change_lines <- geom_vline(xintercept = phase_changes + 0.5, linetype = "dashed")
  }
  
  SCD_graph + phase_change_lines
}


input <- list()
input$behavior <- "Event behavior"
input$freq <- 3
input$dispersion <- 1
input$freq_change <- -50
input$duration <- NA
input$interim_time <- NA
input$duration_change <- NA
input$interim_change <- NA
input$system <- "Frequency counting"
input$interval_length <- 15
input$session_length <- 5
input$design <- "Treatment Reversal"
input$cases <- 3
input$phase_pairs <- 2
input$sessions_TR <- 5
input$sessions_MB <- 30
input$phase_change_list <- "5,10,15"
input$refresh <- NA
input$samples <- 5

design <- phase_design(input$design, input$cases, input$phase_pairs, input$sessions_TR, 
                    input$sessions_MB, input$phase_change_list, input$samples)

dat <- design$dat
dat$Y <- simulate_measurements(dat$trt, input$behavior, 
                               input$freq, input$dispersion, input$freq_change, 
                               input$duration, input$interim_time, input$duration_change, input$interim_change,
                               input$system, input$interval_length, input$session_length)

graph_SCD(dat, input$design, design$phase_changes, input$system)
