library(shiny)

ui <- fluidPage(
  headerPanel('Single-Case Design Simulator'),
  
  fluidRow(

    # Behavioral parameters
    column(4, h3("Behavioral parameters"),
           selectInput("behavior", label = "Behavior class", 
                       choices = c("Event behavior", "State behavior")),
           conditionalPanel(
             condition = "input.behavior=='Event behavior'",
             numericInput("freq", label = "Frequency (per min)", value = NA, min = 0, step = 0.1),
             numericInput("dispersion", label = "Variability", value = 1, min = 0),
             numericInput("freq_change", label = "Percentage change in frequency", value = 0, min = -100, step = 10)
           ),
           conditionalPanel(
             condition = "input.behavior=='State behavior'",
             numericInput("duration", label = "Event duration (seconds)", value = NA, min = 0, step = 1),
             numericInput("interim_time", label = "Interim time (seconds)", value = NA, min = 0, step = 1),
             numericInput("duration_change", label = "Percentage change in event duration", value = 0, min = -100, step = 10),
             numericInput("interim_change", label = "Percentage change in Interim time", value = 0, min = -100, step = 10)
           )
    ),

    # Measurement procedures
    column(4, h3("Measurement procedures"),
           htmlOutput("systemUI"),
           conditionalPanel(
             condition = "input.system=='Momentary time sampling'||input.system=='Partial interval recording'||input.system=='Whole interval recording'",
             numericInput("interval_length", label = "Interval length (seconds)", value = 15, min = 1)
           ),
           numericInput("session_length", label = "Session length (min)", value = 5, min = 1) 
    ),
    
    # Study design
    column(4, h3("Study design"),
           selectInput("design", label = "Study design", choices = c("Treatment Reversal","Multiple Baseline")),
           conditionalPanel(
             condition = "input.design=='Treatment Reversal'",
             numericInput("cases_TR", label = "Number of cases", value = 1, min = 1, max = 12),
             numericInput("phase_pairs", label = "Number of (AB) phases", value = 2, min = 1),
             numericInput("sessions_TR", label = "Sessions per phase", value = 5, min = 1)
           ),
           conditionalPanel(
             condition = "input.design=='Multiple Baseline'",
             numericInput("cases_MB", label = "Number of cases", value = 3, min = 1, max = 12),
             numericInput("sessions_MB", label = "Total number of Sessions", value = 5, min = 1),
             textInput("phase_changes", label = "Phase change times", value = "")
           )
    )
  ),
  
  # Output
  plotOutput('SCDplot'),
  
  actionButton("refresh", label = "Refresh"),
  numericInput("samples", label = "Samples per case", value = 1, min = 1, max = 100)
)

input <- NULL
input$behavior <- "Event behavior"
input$freq <- 3
input$dispersion <- 1
input$freq_change <- -50
input$duration <- NA
input$interim_time <- NA
input$duration_change <- NA
input$interim_change <- NA
input$system <- "Partial interval recording"
input$interval_length <- 15
input$session_length <- 5
input$design <- "Multiple Baseline"
input$cases_TR <- NA
input$phase_pairs <- NA
input$sessions_TR <- NA
input$cases_MB <- 3
input$sessions_MB <- 30
input$phase_changes <- "5, 15, 25"
input$refresh <- NA
input$samples <- 5

phase_design <- function() {
  
}

simulate_measurements <- function() {
  
}

graph_MB <- function() {
  
}

graph_TR <- function() {
  
}

server <- function(input, output) {

  choices <- c("Frequency counting","Continuous recording", 
               "Momentary time sampling","Partial interval recording","Whole interval recording")
  
  output$systemUI <- renderUI( {
    choices_available <- switch(input$behavior,
                                "Event behavior" = choices[c(1,4)],
                                "State behavior" = choices[-1])
    selectInput("system", label = "Measurement system", choices = choices_available)
  })
}

shinyApp(ui = ui, server = server)
