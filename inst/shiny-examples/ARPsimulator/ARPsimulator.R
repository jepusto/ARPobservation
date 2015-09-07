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
    phase <- rep(1:n_alternations, each = n)
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
  impact <- matrix(c(as.numeric(trt[1]==(levels(trt)[-1])), rep(NA, (n - 1) * (nlevels(trt) - 1))),
                   nlevels(trt) - 1, n)
  for (j in 2:n) impact[,j] <- omega * impact[,j - 1] + (1 - omega) * (trt[j]==(levels(trt)[-1]))
  data.frame(t(impact))
}

simulate_measurements <- function(dat, behavior, freq, freq_dispersion, 
                                  duration, interim_time, state_dispersion, 
                                  trt_effect_params,
                                  system, interval_length, session_length) {

  N <- nrow(dat)  
  omega <- 1 - trt_effect_params$immediacy / 100
  group_by(dat, sample, case) %>%
    do(impact(.$trt, omega = omega)) %>%
    ungroup() %>% select(-sample, -case) %>%
    as.matrix() ->
    X_impact
  
  if (behavior == "Event behavior") {
    mu <- rep(0, N)
    lambda <- 1 / (freq * (X_impact %*% trt_effect_params$freq_change / 100 + 1))
    F_event <- F_const()
    F_interim <- F_gam(shape = 1 / freq_dispersion)
  } else {
    mu <- duration * (X_impact %*% trt_effect_params$duration_change / 100 + 1) / 60
    lambda <- interim_time * (X_impact %*% trt_effect_params$interim_change / 100 + 1) / 60
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
    SCD_graph <- SCD_graph + geom_vline(data = phase_dat, aes(xintercept = phase_change), linetype = "dashed")
  } else if (design == "Treatment Reversal") {
    SCD_graph <- SCD_graph + geom_vline(xintercept = phase_changes + 0.5, linetype = "dashed")
  }
  
  SCD_graph
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

