library(shiny)
library(markdown)
library(ARPobservation)
library(dplyr)
library(tidyr)
library(ggplot2)
source("effect_sizes.R")
source("ARPsimulator.R")

trt_effect_UI <- function(k) {
  column(12, 
    conditionalPanel(
      condition = "input.n_trt > 1",
      h4(paste("Treatment",LETTERS[k+1]))
    ),
    conditionalPanel(
       condition = "input.behavior=='Event behavior'",
       numericInput(paste0("freq_change",k), label = "Percentage change in frequency", value = 0, min = -100, step = 10)
    ),
    conditionalPanel(
      condition = "input.behavior=='State behavior'",
      numericInput(paste0("duration_change",k), label = "Percentage change in event duration", value = 0, min = -100, step = 10),
      numericInput(paste0("interim_change",k), label = "Percentage change in Interim time", value = 0, min = -100, step = 10)
    ),
    sliderInput(paste0("immediacy",k), label = "Immediacy of change (%)", 
                min = 0, max = 100, value = 100, step = 5)
  )
}

server <- function(input, output) {

  output$trt_effects_UI <- renderUI({
    lapply(1:input$n_trt, trt_effect_UI)
  })
  
  trt_effect_params <- reactive({
    
    trts <- 1:input$n_trt
    
    freq_change <- rep(0, input$n_trt)
    duration_change <- rep(0, input$n_trt)
    interim_change <- rep(0, input$n_trt)
    immediacy <- rep(0, input$n_trt)
    if (any(grepl("_change",names(input)))) {
      for (t in trts) {
        freq_change[t] <- input[[paste0("freq_change",t)]]
        duration_change[t] <- input[[paste0("duration_change",t)]]
        interim_change[t] <- input[[paste0("interim_change",t)]]
        immediacy[t] <- input[[paste0("immediacy",t)]]
      }
    }

    list(freq_change = freq_change, duration_change = duration_change, 
         interim_change = interim_change, immediacy = immediacy)
  })
    
  choices <- c("Frequency counting","Continuous recording", 
               "Momentary time sampling","Partial interval recording","Whole interval recording")
  
  output$systemUI <- renderUI( {
    choices_available <- switch(input$behavior,
                                "Event behavior" = choices[c(1,4)],
                                "State behavior" = choices[-1])
    selectInput("system", label = "Measurement system", choices = choices_available)
  })
  
  output$cases_UI <- renderUI({
    cases <- 1L + 2 * (input$design == "Multiple Baseline")
    numericInput("cases", label = "Number of cases", value = cases, min = 1)
  })
  
  output$TR_phase_pattern_UI <- renderUI({
    phases <- paste0(LETTERS[rep(1:(input$n_trt + 1), 2)], collapse = "")
    textInput("phase_pattern", label = "Phase pattern", value = phases)
  })
  
  output$MB_phase_change_UI <- renderUI({
    cases <- if (is.null(input$cases)) 1L else input$cases
    n_trt <- if (is.null(input$n_trt)) 1L else input$n_trt
    phase_changes <- trunc(input$sessions_MB * (1:(cases * n_trt)) / (cases * n_trt + 1))
    lapply(1:input$n_trt, function(k) {
      lab <- if (input$n_trt > 1) paste0("Phase change times (Trt ", LETTERS[k+1],")") else "Phase change times"
      textInput(paste0("phase_change_list",k), 
                label = lab, 
                value = paste(phase_changes[cases * (k - 1) + 1:cases], collapse = ", "))  
    })
  })
  
  MB_phase_changes <- reactive({
    phase_change_list <- c()
    for (i in 1:(input$n_trt)) {
      phase_change_list[i] <- input[[paste0("phase_change_list",i)]]
    }
    phase_change_list
  })
  
  sim_dat <- eventReactive(c(input$outputPanel, input$simulateGraph, input$simulateES), {
    cases <- if (is.null(input$cases)) 1L else input$cases
    system <- if (is.null(input$system)) {
      if (input$behavior=="Event behavior") choices[1] else choices[2]
    } else {
      input$system
    }
    phase_pattern <- if (is.null(input$phase_pattern)) "ABAB" else input$phase_pattern

    phase_changes <- get_phase_changes(input$design, input$sessions_TR, phase_pattern, 
                                       MB_phase_changes(), input$cases)
    
    samples <- ifelse(input$outputPanel == "SCD Graph", input$samplesGraph, input$samplesES)
    
    dat <- phase_design(input$design, input$n_trt, cases, phase_pattern, 
                        input$sessions_TR, input$sessions_MB, phase_changes, 
                        input$n_alternations, input$randomize_AT, samples)
    
    dat <- simulate_measurements(dat, input$behavior, input$freq, input$freq_dispersion, 
                                 input$duration, input$interim_time, input$state_dispersion,
                                 trt_effect_params(), 
                                 system, input$interval_length, input$session_length)
    
    height <- max(300, 150 * cases)
    list(dat = dat, design = input$design, phase_changes = phase_changes, 
         system = system, height_SCD = height)
  })
  
  output$SCDplot <- renderPlot({
    if (input$simulateGraph > 0 | input$simulateES > 0) {
      with(sim_dat(), graph_SCD(dat, design, phase_changes, system, input$showtruth))
    }
  }, height = function() sim_dat()$height_SCD)

  output$phase_pre_UI <- renderUI({
    selectInput("phase_pre", label = "Pre phase", choices = LETTERS[1:(input$n_trt + 1)])
  })
  
  output$phase_post_UI <- renderUI({
    selectInput("phase_post", label = "Post phase", 
                choices = setdiff(LETTERS[1:(input$n_trt + 1)], input$phase_pre))
  })
  

  ES_dat <- reactive({
    if (input$simulateGraph > 0 | input$simulateES > 0) {
      ES_dat <- calculate_ES(sim_dat()$dat, input$phase_pre, input$phase_post, 
                             input$effect_size, input$improvement)
    } else {
      ES_dat <- NULL
    }
    ES_dat
  })
  
  output$ESplot <- renderPlot({
    if (input$simulateGraph > 0 | input$simulateES > 0) {
      graph_ES(ES_dat(), input$effect_size, input$showAvgES)
    }
    }, height = 400)
  
}
