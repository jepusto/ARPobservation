library(shiny)
library(ARPobservation)
library(ggplot2)
source("ARPsimulator.R")

ui <- fluidPage(
  headerPanel('Alternating Renewal Process Simulator'),
  
  fluidRow(

    # Behavioral parameters
    column(4, h3("Behavioral parameters"),
           selectInput("behavior", label = "Behavior class", 
                       choices = c("Event behavior", "State behavior")),
           conditionalPanel(
             condition = "input.behavior=='Event behavior'",
             numericInput("freq", label = "Frequency (per min)", value = NA, min = 0, step = 0.1),
             numericInput("dispersion", label = "Variability", value = 1, min = 0.05, step = 0.05),
             numericInput("freq_change", label = "Percentage change in frequency", value = 0, min = -100, step = 10)
           ),
           conditionalPanel(
             condition = "input.behavior=='State behavior'",
             numericInput("duration", label = "Event duration (seconds)", value = NA, min = 0, step = 1),
             numericInput("interim_time", label = "Interim time (seconds)", value = NA, min = 0, step = 1),
             numericInput("duration_change", label = "Percentage change in event duration", value = 0, min = -100, step = 10),
             numericInput("interim_change", label = "Percentage change in Interim time", value = 0, min = -100, step = 10)
           ),
           sliderInput("immediacy", label = "Immediacy of change (%)", 
                       min = 0, max = 100, value = 100, step = 1)
    ),

    # Study design
    column(4, h3("Study design"),
           selectInput("design", label = "Study design", choices = c("Treatment Reversal","Multiple Baseline")),
           htmlOutput("cases_UI"),
           conditionalPanel(
             condition = "input.design=='Treatment Reversal'",
             numericInput("phase_pairs", label = "Number of (AB) phases", value = 2, min = 1),
             numericInput("sessions_TR", label = "Sessions per phase", value = 5, min = 1)
           ),
           conditionalPanel(
             condition = "input.design=='Multiple Baseline'",
             numericInput("sessions_MB", label = "Total number of Sessions", value = 20, min = 1),
             htmlOutput("MB_phase_change_UI")
           )
    ),
    
    # Measurement procedures
    column(4, h3("Measurement procedures"),
           numericInput("session_length", label = "Session length (min)", value = 10, min = 1),
           htmlOutput("systemUI"),
           conditionalPanel(
             condition = "input.system=='Momentary time sampling'||input.system=='Partial interval recording'||input.system=='Whole interval recording'",
             numericInput("interval_length", label = "Interval length (seconds)", value = 15, min = 1)
           )
    )
  ),

  fluidRow(column(12, 
                  br(),
                  hr(),
                  br())),
  
  tabsetPanel(id = "outputPanel", type = "tabs",
              tabPanel("Graph",  
                
                fluidRow(column(12, br())),
                plotOutput('SCDplot', height = "auto"),
  
                fluidRow(
                  column(3, numericInput("samples", label = "Samples per case", value = 1, min = 1, max = 100)),
                  column(6, align = "center", actionButton("simulate", label = "Simulate!")),
                  column(3, checkboxInput("showtruth", label = "Show true trend lines", value = FALSE))
                )
              ),
              tabPanel("Effect sizes"))
)

server <- function(input, output) {

  choices <- c("Frequency counting","Continuous recording", 
               "Momentary time sampling","Partial interval recording","Whole interval recording")
  
  output$systemUI <- renderUI( {
    choices_available <- switch(input$behavior,
                                "Event behavior" = choices[c(1,4)],
                                "State behavior" = choices[-1])
    selectInput("system", label = "Measurement system", choices = choices_available)
  })

  output$cases_UI <- renderUI({
    cases <- 1 + 2 * (input$design == "Multiple Baseline")
    numericInput("cases", label = "Number of cases", value = cases, min = 1)
  })
  
  output$MB_phase_change_UI <- renderUI({
      cases <- if (is.null(input$cases)) 1 else input$cases
      phase_changes <- trunc(input$sessions_MB * (1:cases) / (cases + 1))
      phase_change_list <- paste(phase_changes, collapse = ", ")
      textInput("phase_change_list", label = "Phase change times", value = phase_change_list)  
  })
  
  sim_dat <- eventReactive(input$simulate, {
    phase_changes <- get_phase_changes(input$design, input$sessions_TR, input$phase_pairs, 
                                       input$phase_change_list, input$cases)
    dat <- phase_design(input$design, input$cases, input$phase_pairs, input$sessions_TR, 
                           input$sessions_MB, phase_changes, input$samples)
    dat <- simulate_measurements(dat, input$behavior, 
                                   input$freq, input$dispersion, input$freq_change, 
                                   input$duration, input$interim_time, input$duration_change,
                                   input$interim_change, input$immediacy, 
                                   input$system, input$interval_length, input$session_length)
    height <- max(400, 150 * input$cases)
    list(dat = dat, design = input$design, phase_changes = phase_changes, system = input$system, height = height)
  })

  output$SCDplot <- renderPlot({
    with(sim_dat(), graph_SCD(dat, design, phase_changes, system, input$showtruth))
  }, height = function() sim_dat()$height)
  
}

shinyApp(ui = ui, server = server)
