library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)

# Output function for DataTables
dt_output <- function(title, id) {
  box(
    title = title,
    status = "primary",
    solidHeader = TRUE,
    width = 12,
    DTOutput(id)
  )
}

# Render function for DataTables
render_dt <- function(data, editable = 'cell', server = TRUE, ...) {
  renderDT(data, selection = 'none', server = server, editable = editable, ...)
}

# Define the origin host directory and stages
origin_host_dir <- '../../.'
xStage <- c("I", "II", "III", "IV")

# Read the data
dwell_model_group_df <- read_tsv(sprintf("%s/data/20200728_dwell_time_groups.tsv", origin_host_dir))
dwell_model_timing_df <- read_tsv(sprintf("%s/data/20200728_dwell_group_timing.tsv", origin_host_dir))

# UI function
ui <- dashboardPage(
  dashboardHeader(title = "Dwell Time Model"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidPage(
      box(
        title = "Introduction",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p("This page will guide you through the process of generating cancer dwell time for each cancer and stage.")
      ),
      box(
        title = "Table 1: Cancer Group",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p("For Table 1, double click to edit", strong('group'), ", and click anywhere to save", style = "color:gray"),
        strong("You can only edit the group column"),
        dt_output('Table 1: Cancer group', 'cancer_group')
      ),
      box(
        title = "Table 2: Dwell Time",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p("For Table 2, double click to edit", strong('dwell'), ", and click anywhere to save", style = "color:gray"),
        strong("You can only edit the dwell column"),
        dt_output('Table 2: Dwell time', 'dwell_time')
      ),
      box(
        title = "Scenarios",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p('Scenarios | 0: Maximum interception scenario (MIS) 1: Very slow, 2: Slow, 3: Fast, 4: Aggressively fast')
      ),
      box(
        title = "Table 3: All Data",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        downloadButton("downloadTable3", "Download Table 3"),
        dataTableOutput("dwell_model_all_df")
      ),
      box(
        title = "Dictionary of Variables",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p(HTML('scenario: Unique id for dwell time scenario <br>
               Stage: Stage to which this dwell time applies <br>
               group: dwell group to which a cancer belongs <br>
               dwell: Mean time in years under this scenario <br>
               Cancer: Cancer sensitivity group <br>
               dwell_group: dwell group assigned to this cancer type'))
      )
    )
  )
)

# Server function
server <- function(input, output) {
  
  dwell_model_all_df <- function() {
    updated_df <- dwell_model_group_df %>%
      rename(dwell_group = group) %>%
      full_join(dwell_model_timing_df %>% rename(dwell_group = group)) %>%
      mutate(number_stage = match(Stage, xStage))
    return(updated_df)
  }
  
  output$cancer_group <- render_dt(dwell_model_group_df, list(target = 'cell', server = TRUE), options = list(pageLength = 10))
  
  observeEvent(input$cancer_group_cell_edit, {
    dwell_model_group_df <<- editData(dwell_model_group_df, input$cancer_group_cell_edit, 'cancer_group')
    output$dwell_model_all_df <- renderDataTable(dwell_model_all_df(), options = list(pageLength = 10))
  })
  
  output$dwell_time <- render_dt(dwell_model_timing_df, list(target = 'cell', disable = list(columns = c(1, 2, 3)), server = TRUE), options = list(pageLength = 10))
  
  observeEvent(input$dwell_time_cell_edit, {
    dwell_model_timing_df <<- editData(dwell_model_timing_df, input$dwell_time_cell_edit, 'dwell_time')
    output$dwell_model_all_df <- renderDataTable(dwell_model_all_df(), options = list(pageLength = 10))
  })
  
  output$dwell_model_all_df <- renderDataTable(dwell_model_all_df(), options = list(pageLength = 10))
  
  # Download final table
  output$downloadTable3 <- downloadHandler(
    filename = function() {
      "dwell_model_all_df.csv"
    },
    content = function(file) {
      write.csv(dwell_model_all_df(), file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui, server)