library(shiny)

ui <- fluidPage(
  headerPanel('Single-Case Design Simulator'),
  
  fluidRow(

    # Behavioral parameters
    column(4, h3("Behavioral parameters"),
           selectInput("behavior", label = "Behavior class", 
                       choices = c("Event behavior", "State behavior")),
           numericInput("freq", label = "Frequency (per min)", value = NA, min = 0, step = 0.1),
           numericInput("dispersion", label = "Variability", value = 1, min = 0),
           numericInput("freq_change", label = "Percentage change in frequency", value = 0, min = -100)
    ),

    # Measurement procedures
    column(4, h3("Measurement procedures"),
           selectInput("system", label = "Measurement system", 
                       choices = c("Frequency counting","Partial interval recording","Continuous recording", 
                                   "Momentary time sampling", "Whole interval recording")),
           numericInput("session_length", label = "Session length (min)", value = 5, min = 1) 
    ),
    
    # Study design
    column(4, h3("Study design"),
           selectInput("design", label = "Study design", choices = c("Treatment Reversal","Multiple Baseline")),
           conditionalPanel(
             condition = "input.design=='Treatment Reversal'",
             numericInput("cases", label = "Number of cases", value = 1, min = 1, max = 12),
             numericInput("phase_pairs", label = "Number of (AB) phases", value = 2, min = 1),
             numericInput("sessions", label = "Sessions per phase", value = 5, min = 1)
           ),
           numericInput("samples", label = "Samples per case", value = 1, min = 1, max = 100)
    )
  ),
  
  # Output
  plotOutput('SCDplot'),
  actionButton("refresh", label = "Refresh")
)

server <- function(input, output) {}

shinyApp(ui = ui, server = server)
