library(ARPobservation)
library(dplyr)
library(ggplot2)
rm(list=ls())
source("inst/shiny-examples/ARPsimulator/effect_sizes.R")
source("inst/shiny-examples/ARPsimulator/ARPsimulator.R")

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
input$session_length <- 10

# Study design
input$design <- "Treatment Reversal"
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
input$samples <- 2
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

dat <- phase_design(input$design, input$n_trt, input$cases, input$phase_pattern, 
                    input$sessions_TR, input$sessions_MB, phase_changes, 
                    input$n_alternations, input$randomize_AT, input$samples)

dat <- simulate_measurements(dat, input$behavior, input$freq, input$freq_dispersion, 
                             input$duration, input$interim_time, input$state_dispersion,
                             trt_effect_params, 
                             input$system, input$interval_length, input$session_length)

sim_dat <- list(dat = dat, design = input$design, phase_changes = phase_changes, 
                system = input$system)

with(sim_dat, graph_SCD(dat, design, phase_changes, system, input$showtruth))

with(sim_dat, graph_ES(dat, input$effect_size, input$improvement))
