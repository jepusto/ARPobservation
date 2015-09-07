#------------------------------
# Design matrix
#------------------------------

get_phase_changes <- function(design, sessions_TR, phase_pattern, MB_phase_changes, cases) {
  if (design=="Treatment Reversal") {
    phase_labels <- strsplit(phase_pattern, "")[[1]]
    phase_changes <- sessions_TR * (1:(length(phase_labels) - 1))
  } else {
    phase_changes <- lapply(MB_phase_changes, function(x) as.numeric(strsplit(x, split = ",")[[1]]))
    phase_changes <- sapply(phase_changes, function(x) rep(x, length.out = cases))
    phase_changes <- as.list(data.frame(t(phase_changes)))
    names(phase_changes) <- NULL
  }
  return(phase_changes)
}

phase_design <- function(design, n_trt, cases, phase_pattern, sessions_TR, 
                         sessions_MB, phase_changes, n_alternations, randomize_AT, samples) {
  
  case_names <- paste("Case", LETTERS[1:cases])
  
  sample_nr <- 1:samples

  if (design=="Treatment Reversal") {
    phase_labels <- strsplit(phase_pattern, "")[[1]]
    session_nr <- 1:(length(phase_labels) * sessions_TR)
    phase <- rep(1:length(phase_labels), each = sessions_TR)
    trt <- rep(phase_labels, each = sessions_TR)
  } else if (design=="Multiple Baseline") {
    session_nr <- 1:sessions_MB
    phase <- unlist(lapply(phase_changes, function(t) 
      rep(1:(n_trt + 1), times = diff(c(0,t,sessions_MB)))))
    trt <- LETTERS[phase]
  } else {
    n <- n_trt + 1
    session_nr <- 1:(n * n_alternations)
    if (randomize_AT) {
      trt <- as.vector(replicate(n_alternations, sample(LETTERS[1:n])))
    } else {
      trt <- rep(LETTERS[1:(n_trt + 1)], n_alternations)
    }
  }
  
  dat <- expand.grid(session = session_nr, case = case_names, sample = sample_nr)
  dat$phase <- factor(phase)
  dat$trt <- factor(trt)
  return(dat)
}

#----------------------------
# Generate outcome data
#----------------------------

impact <- function(trt, omega) {
  n <- length(trt)
  impact <- c(as.numeric(trt[1]=="Treat"), rep(NA, n - 1))
  for (j in 2:n) impact[j] <- omega * impact[j - 1] + (1 - omega) * (trt[j]=="Treat")
  impact
}

simulate_measurements <- function(dat, behavior, freq, freq_dispersion, 
                                  duration, interim_time, state_dispersion, 
                                  freq_change, duration_change, interim_change, immediacy,
                                  system, interval_length, session_length) {

  N <- nrow(dat)  
  omega <- 1 - immediacy / 100
  impact <- with(dat, unlist(tapply(trt, list(case, sample), impact, omega = omega)))
  
  if (behavior == "Event behavior") {
    mu <- rep(0, N)
    lambda <- 1 / (freq * (impact * freq_change / 100 + 1))
    F_event <- F_const()
    F_interim <- F_gam(shape = 1 / freq_dispersion)
  } else {
    mu <- duration * (impact * duration_change / 100 + 1) / 60
    lambda <- interim_time * (impact * interim_change / 100 + 1) / 60
    F_event <- F_gam(shape = 1 / state_dispersion)
    F_interim <- F_gam(shape = 1 / state_dispersion)
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

library(ARPobservation)
library(dplyr)
library(ggplot2)
source("inst/shiny-examples/ARPsimulator/effect_sizes.R")

input <- list()

# baseline behavior 
input$behavior <- "Event behavior"
input$freq <- 3
input$freq_dispersion <- 1
input$duration <- 30
input$interim_time <- 90
input$state_dispersion <- 1

# behavior change
input$n_trt <- 3
input$freq_change1 <- -50
input$duration_change1 <- NA
input$interim_change1 <- NA
input$immediacy1 <- 20
input$freq_change2 <- -80
input$duration_change2 <- NA
input$interim_change2 <- NA
input$immediacy2 <- 20
input$freq_change3 <- -90
input$duration_change3 <- NA
input$interim_change3 <- NA
input$immediacy3 <- 20

# Measurement system
input$system <- "Frequency counting"
input$interval_length <- 15
input$session_length <- 100

# Study design
input$design <- "Treatment reversal"
input$cases <- 3
input$phase_pattern <- "ABCDABCD"
input$sessions_TR <- 5
input$sessions_MB <- 30
input$phase_change_list1 <- "5,10,15"
input$phase_change_list2 <- "12,17,22"
input$phase_change_list3 <- "19,24,28"
input$n_alternations <- 3
input$randomize_AT <- TRUE

# Miscellaneous
input$refresh <- NA
input$samples <- 1
input$showtruth <- TRUE
input$effect_size <- "NAP"
input$improvement <- 2

# Treatment effect parameters (reactive)
trts <- 1:input$n_trt
freq_change <- unlist(input[paste0("freq_change",trts)])
duration_change <- unlist(input[paste0("duration_change",trts)])
interim_change <- unlist(input[paste0("interim_change",trts)])
immediacy <- unlist(input[paste0("immediacy",trts)])

trt_effect_params <- list(freq_change = freq_change, duration_change = duration_change, 
                          interim_change = interim_change, immediacy = immediacy)


# MB phase changes (reactive)
MB_phase_changes <- unlist(input[paste0("phase_change_list",1:input$n_trt)])


phase_changes <- get_phase_changes(input$design, input$sessions_TR, input$phase_pattern, 
                                   MB_phase_changes, input$cases)

dat <- phase_design(input$design, input$n_trt, input$cases, input$phase_pairs, input$sessions_TR, 
                    input$sessions_MB, phase_changes, 
                    input$n_alternations, input$randomize_AT, input$samples)

dat <- simulate_measurements(dat, input$behavior, 
                               input$freq, input$dispersion, input$freq_change, 
                               input$duration, input$interim_time, input$duration_change, 
                               input$interim_change, input$immediacy, 
                               input$system, input$interval_length, input$session_length)

sim_dat <- list(dat = dat, design = input$design, phase_changes = phase_changes, 
                system = input$system)

with(sim_dat, graph_SCD(dat, design, phase_changes, system, input$showtruth))

with(sim_dat, graph_ES(dat, input$effect_size, input$improvement))
