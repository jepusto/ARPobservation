library(shiny)

ui <- navbarPage(title = "Alternating Renewal Process Simulator",
   tabPanel("Simulator",
            
      fluidRow(
        
        # Baseline behavioral parameters
        column(3, h3("Baseline behavior"),
               selectInput("behavior", label = "Behavior class", 
                           choices = c("Event behavior", "State behavior")),
               conditionalPanel(
                 condition = "input.behavior=='Event behavior'",
                 numericInput("freq", label = "Frequency (per min)", value = 1, min = 0, step = 0.1),
                 numericInput("freq_dispersion", label = "Variability", value = 1, min = 0.05, step = 0.05)
               ),
               conditionalPanel(
                 condition = "input.behavior=='State behavior'",
                 numericInput("duration", label = "Event duration (seconds)", value = 30, min = 0, step = 1),
                 numericInput("interim_time", label = "Interim time (seconds)", value = 60, min = 0, step = 1),
                 numericInput("state_dispersion", label = "Variability", value = 1, min = 0.05, step = 0.05)
               )
        ),
        
        # Behavior change parameters
        column(3, h3("Behavior change"),
               numericInput("n_trt", label = "Number of treatments", value = 1, min = 1, max = 4),
               htmlOutput("trt_effects_UI")
        ),
        
        # Measurement procedures
        column(3, h3("Measurement procedures"),
               numericInput("session_length", label = "Session length (min)", value = 10, min = 1),
               htmlOutput("systemUI"),
               conditionalPanel(
                 condition = "input.system=='Momentary time sampling'||input.system=='Partial interval recording'||input.system=='Whole interval recording'",
                 numericInput("interval_length", label = "Interval length (seconds)", value = 15, min = 1)
               )
        ),
        
        # Study design
        column(3, h3("Study design"),
               selectInput("design", label = "Study design", choices = c("Treatment Reversal","Multiple Baseline","Alternating Treatment")),
               htmlOutput("cases_UI"),
               conditionalPanel(
                 condition = "input.design=='Treatment Reversal'",
                 htmlOutput("TR_phase_pattern_UI"),
                 numericInput("sessions_TR", label = "Sessions per phase", value = 5, min = 1)
               ),
               conditionalPanel(
                 condition = "input.design=='Multiple Baseline'",
                 numericInput("sessions_MB", label = "Total number of sessions", value = 20, min = 1),
                 htmlOutput("MB_phase_change_UI")
               ),
               conditionalPanel(
                 condition = "input.design=='Alternating Treatment'",
                 numericInput("n_alternations", label = "Number of alternations", value = 5, min = 1),
                 checkboxInput("randomize_AT", label = "Randomize treatment order", value = TRUE)
               )
        )
      ),
      
      fluidRow(column(12, 
                      hr(),
                      h3("Results"))),
      
      tabsetPanel(id = "outputPanel", type = "tabs",
                  tabPanel("SCD Graph",  
                           column(12, br()),
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          numericInput("samplesGraph", label = "Samples per case", value = 1, min = 1, max = 100),
                                          checkboxInput("showtruth", label = "Show true trend lines", value = FALSE),
                                          column(12, align = "center", actionButton("simulateGraph", label = "Simulate!")),
                                          br()
                             ),
                             mainPanel(width = 9,
                                       plotOutput('SCDplot', height = "auto")
                             )
                           )
                  ),
                  tabPanel("Effect sizes",
                           column(12, br()),
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          conditionalPanel(
                                            condition = "input.n_trt > 1",
                                            column(6, htmlOutput("phase_pre_UI")),
                                            column(6, htmlOutput("phase_post_UI"))
                                          ),
                                          selectInput("effect_size", label = "Effect size measure", 
                                                      choices = c("PND","PEM","PAND","IRD","NAP","Tau","Within-case SMD")),
                                          conditionalPanel(
                                            condition = "input.effect_size != 'Within-case SMD'",
                                            radioButtons("improvement", label = "Direction of improvement", choices = list("increase" = 1, "decrease" = 2), selected = 1)
                                          ),
                                          numericInput("samplesES", label = "Samples per case", value = 100, min = 1, max = 1000),
                                          checkboxInput("showAvgES", label = "Show average", value = FALSE),
                                          column(12, align = "center", actionButton("simulateES", label = "Simulate!")),
                                          br()
                             ),
                             mainPanel(width = 9,
                                       plotOutput('ESplot', height = "auto")
                             )
                           )
                  )
      )
   ),
   tabPanel("Help",
      navlistPanel(widths = c(3,9),
         tabPanel("Overview", includeMarkdown("markdown/Overview.md")),
         tabPanel("Baseline behavior", includeMarkdown("markdown/Behavioral_parameters.md")),
         tabPanel("Behavior change", includeMarkdown("markdown/Behavior_change.md")),
         tabPanel("Measurement procedures", includeMarkdown("markdown/Measurement_procedures.md")),
         tabPanel("Study design features", includeMarkdown("markdown/Study_design.md")),
         tabPanel("Single-case graph", includeMarkdown("markdown/SCD_graph.md")),
         tabPanel("Effect size graph", includeMarkdown("markdown/ES_graph.md"))
      )
   ),
   tabPanel("About",
      navlistPanel(widths = c(3,9),
         tabPanel("ARPsimulator", includeMarkdown("markdown/ARPsimulator.md")),
         tabPanel("Accessing the simulator", includeMarkdown("markdown/Accessing_ARPsimulator.md"))
      )
   )
)
