#------------------------------
# Design matrix
#------------------------------

get_phase_changes <- function(design, sessions_TR, phase_pairs, phase_change_list, cases) {
  if (design=="Treatment Reversal") {
    phase_changes <- sessions_TR * (1:(2 * phase_pairs - 1))
  } else {
    phase_changes <- as.numeric(strsplit(phase_change_list, split = ",")[[1]])
    phase_changes <- rep(phase_changes, length.out = cases)
  }
  return(phase_changes)
}

phase_design <- function(design, cases, phase_pairs, sessions_TR, 
                         sessions_MB, phase_changes, samples) {
  
  case_names <- paste("Case", LETTERS[1:cases])
  
  if (design=="Treatment Reversal") {
    sample_nr <- 1:samples
    session_nr <- 1:(phase_pairs * 2 * sessions_TR)
    phase <- rep(1:(2 * phase_pairs), each = sessions_TR)
    trt <- rep(rep(c("Base","Treat"), each = sessions_TR), times = phase_pairs)
  } else {
    sample_nr <- 1:samples
    session_nr <- 1:sessions_MB
    phase <- unlist(lapply(phase_changes, function(t) rep(1:2, times = c(t, sessions_MB - t))))
    trt <- c("Base","Treat")[phase]
  }
  
  dat <- expand.grid(session = session_nr, case = case_names, sample = sample_nr)
  dat$phase <- factor(phase)
  dat$trt <- factor(trt)
  return(dat)
}

#----------------------------
# Generate outcome data
#----------------------------

go <- function(behavior, freq, duration, interim_time) {
  (behavior == "Event behavior" & !is.na(freq)) | 
       (behavior == "State behavior" & !is.na(duration) & !is.na(interim_time))
}

impact <- function(trt, omega) {
  n <- length(trt)
  impact <- c(as.numeric(trt[1]=="Treat"), rep(NA, n - 1))
  for (j in 2:n) impact[j] <- omega * impact[j - 1] + (1 - omega) * (trt[j]=="Treat")
  impact
}

simulate_measurements <- function(dat, behavior, freq, dispersion, freq_change, 
                                  duration, interim_time, duration_change, 
                                  interim_change, immediacy,
                                  system, interval_length, session_length) {

  N <- nrow(dat)  
  omega <- 1 - immediacy / 100
  impact <- with(dat, unlist(tapply(trt, list(case, sample), impact, omega = omega)))
  
  if (behavior == "Event behavior") {
    mu <- rep(0, N)
    lambda <- 1 / (freq * (impact * freq_change / 100 + 1))
    F_event <- F_const()
    F_interim <- F_gam(shape = 1 / dispersion)
  } else {
    mu <- duration * (impact * duration_change / 100 + 1) / 60
    lambda <- interim_time * (impact * interim_change / 100 + 1) / 60
    F_event <- F_exp()
    F_interim <- F_exp()
  }
  
  BS <- r_behavior_stream(n = N, mu = mu, lambda = lambda,
                          F_event = F_event, F_interim = F_interim, 
                          stream_length = session_length)
  
  # compute true trend lines
  if (system == "Frequency counting") {
    dat$truth <- session_length / lambda
  } else {
    dat$truth <- 100 * mu / (mu + lambda)
  }
  
  # compute observed data
  dat$Y <- switch(system,
              "Frequency counting" = event_counting(BS),
              "Continuous recording" = 100 * continuous_duration_recording(BS),
              "Momentary time sampling" = 100 * momentary_time_recording(BS, interval_length / 60),
              "Partial interval recording" = 100 * interval_recording(BS, interval_length / 60, partial = TRUE),
              "Whole interval recording" = 100 * interval_recording(BS, interval_length / 60, partial = FALSE)
  )
  return(dat)
}


#----------------------------
# Create SCD graph
#----------------------------

graph_SCD <- function(dat, design, phase_changes, system, showtruth) {
  
  Y_lab <- switch(system,
                  "Frequency counting" = "Events",
                  "Continuous recording" = "% Session time",
                  "Momentary time sampling" = "% Moments",
                  "Partial interval recording" = "% Partial intervals",
                  "Whole interval recording" = "% Whole intervals"
  )
  Y_range <- c(0, ifelse(system=="Frequency counting", max(dat$Y) + 1, 100))
  
  samples <- max(dat$sample)
  
  SCD_graph <- ggplot(dat, aes(session, Y, color = trt, group = interaction(phase, sample))) + 
    geom_point(alpha = samples^(-1/4)) + 
    geom_line(alpha = samples^(-1/2)) + 
    coord_cartesian(ylim = Y_range) + 
    theme_bw() + theme(legend.position = "bottom") + 
    labs(color = "Phase", y = Y_lab) 
  
  if (nlevels(dat$case) > 1) {
    SCD_graph <- SCD_graph + facet_grid(case ~ .)
  }

  if (showtruth) {
    SCD_graph <- SCD_graph + geom_line(data = subset(dat, sample==1), aes(session, truth), size = 1.0)
  }
  
  if (design == "Multiple Baseline") {
    phase_dat <- data.frame(case = levels(dat$case), phase_change = phase_changes + 0.5)
    phase_change_lines <- geom_vline(data = phase_dat, aes(xintercept = phase_change), linetype = "dashed")
  } else {
    phase_change_lines <- geom_vline(xintercept = phase_changes + 0.5, linetype = "dashed")
  }
  
  SCD_graph + phase_change_lines
}

#----------------------------
# Create effect size graph
#----------------------------

ES_choices <- c("PND","PEM","PAND","IRD","NAP","Tau","Within-case SMD")

graph_ES <- function(dat, effect_size, improvement, showAvgES) {

  ES_function <- switch(effect_size,
                        PND = PND,
                        PEM = PEM,
                        PAND = PAND,
                        IRD = IRD,
                        NAP = NAP,
                        Tau = Tau,
                        "Within-case SMD" = SMD)

  # calculate effect sizes
  group_by(dat, case, sample) %>%
    summarize(ES = ES_function(data = Y, phase = trt, base_phase = "Base", increase = improvement==1)) %>%
    ungroup() -> ES_dat
  
  # graph effect sizes
  X_range <- switch(effect_size,
                    PND = c(0,100),
                    PEM = c(0,100),
                    PAND = c(0,100),
                    IRD = c(-1,1),
                    NAP = c(0,100),
                    Tau = c(-1,1),
                    "Within-case SMD" = range(ES_dat$ES))
  cases <- nlevels(ES_dat$case)
  legend_position <- if (cases > 1) "bottom" else "none"
  ES_graph <- ggplot(ES_dat, aes(ES, fill = case)) + 
    geom_density(alpha = 1 / max(2,cases)) + 
    coord_cartesian(xlim = X_range) + 
    theme_minimal() + theme(legend.position = legend_position) + 
    labs(x = effect_size, y = "distribution", fill = "") 
  
  if (showAvgES) {
    ES_avg <- group_by(ES_dat, case) %>% summarize(ESavg = mean(ES))
    ES_graph <- ES_graph + 
      geom_vline(data = ES_avg, aes(xintercept = ESavg, color = case), size = 1, linetype = "dashed")
  } 
  
  ES_graph
} 

# library(ARPobservation)
# library(dplyr)
# library(ggplot2)
# source("inst/shiny-examples/ARPsimulator/effect sizes.R")
# 
# input <- list()
# input$behavior <- "Event behavior"
# input$freq <- 3
# input$dispersion <- 1
# input$freq_change <- -50
# input$duration <- NA
# input$interim_time <- NA
# input$duration_change <- NA
# input$interim_change <- NA
# input$immediacy <- 20
# input$system <- "Frequency counting"
# input$interval_length <- 15
# input$session_length <- 100
# input$design <- "Multiple Baseline"
# input$cases <- 3
# input$phase_pairs <- 2
# input$sessions_TR <- 5
# input$sessions_MB <- 30
# input$phase_change_list <- "5,10,15"
# input$refresh <- NA
# input$samples <- 50
# input$showtruth <- TRUE
# input$effect_size <- "NAP"
# input$improvement <- 2
# 
# phase_changes <- get_phase_changes(input$design, input$sessions_TR, input$phase_pairs, 
#                                    input$phase_change_list, input$cases)
# 
# dat <- phase_design(input$design, input$cases, input$phase_pairs, input$sessions_TR, 
#                     input$sessions_MB, phase_changes, input$samples)
# 
# dat <- simulate_measurements(dat, input$behavior, 
#                                input$freq, input$dispersion, input$freq_change, 
#                                input$duration, input$interim_time, input$duration_change, 
#                                input$interim_change, input$immediacy, 
#                                input$system, input$interval_length, input$session_length)
# 
# sim_dat <- list(dat = dat, design = input$design, phase_changes = phase_changes, 
#                 system = input$system)
# 
# with(sim_dat, graph_SCD(dat, design, phase_changes, system, input$showtruth))
# 
# with(sim_dat, graph_ES(dat, input$effect_size, input$improvement))
