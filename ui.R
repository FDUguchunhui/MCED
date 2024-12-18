library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(ggalluvial)
library(ggplot2)
library(patchwork)

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



# ui.R ----------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Multi-cancer early detection Interception Model", titleWidth = 450),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction", icon = icon("info-circle")),
      menuItem("Dwell Time", tabName = "dwell_time", icon = icon("clock")),
      menuItem("Cancer Group", tabName = "cancer_group", icon = icon("group")),
      menuItem("Scenarios", tabName = "scenarios", icon = icon("list")),
      menuItem("Sensitivity Table", tabName = "sensitivity_table", icon = icon("table")),
      menuItem("Model Slip Parameters", tabName = "model_slip_parameters", icon = icon("sliders")),
      menuItem("Simulate Performance", tabName = "simulate_performance", icon = icon("chart-line")),
      menuItem("Flow Diagrams", tabName = "flow_diagrams", icon = icon("project-diagram")),
      menuItem("Utility Sankey Plot", tabName = "utility_sankey_plot", icon = icon("sitemap")),
      menuItem("Sensitivity Analysis", tabName = "sensitivity_analysis", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "introduction",
              fluidPage(
                box(
                  title = "Introduction",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  p("This page will guide you through the process of simulated intercepted cancer cases and deaths at each stage and cancer type.")
                )
              )
      ),
      tabItem(tabName = "dwell_time",
              fluidPage(
                box(
                  title = "Step 1: Define different dwell group and dwell time",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  p("Double click to edit the table", strong('dwell'), ", and click anywhere to save", style = "color:gray"),
                  strong("You can only edit the dwell column"),
                  dt_output('Dwell time', 'dwell_time')
                )
              )
      ),
      tabItem(tabName = "cancer_group",
              fluidPage(
                box(
                  title = "Step 2: select dwell group for cancer",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  p("Double click to edit the table", strong('group'), ", and click anywhere to save", style = "color:gray"),
                  strong("You can only edit the group column"),
                  dt_output('Cancer group', 'cancer_group')
                )
              )
      ),
      tabItem(tabName = "scenarios",
              fluidPage(
                box(
                  title = "Scenarios",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  p('Scenarios | 0: Maximum interception scenario (MIS) 1: Very slow, 2: Slow, 3: Fast, 4: Aggressively fast')
                ),
                box(
                  title = "Final dwell time table",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  downloadButton("downloadDwellTable", "Download Table"),
                  dataTableOutput("dwell_model_all_df"),
                  p(HTML('The dwell time for each cancer is defined by which dwell group it belong to and the progress scenario
                   this table is to show the final configuration of the dwell time for each cancer type. <br>')),
                  p(HTML('Cancer: Cancer sensitivity group <br>
                   dwell_group: the pre-defined dwell group, each dwell group has a fixed dwell time for stage 1-4 <br>
                   scenario: 1."VSlow": very slow, 2. slow 3. Fast 4. "AggFast": aggressively fast  <br>
                   stage: the character stage name <br>
                   number_stage: the numerical representative of the stage I: 1, II: 2, III: 3, IV: 4 <br>'))
                )
              )
      ),
      tabItem(tabName = "sensitivity_table",
              fluidPage(
                box(
                  title = "Step 3: Isotonic Regression Adjusted Sensitivity Table",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  p('Assumption: Sensitivity estimates are nondecreasing by increasing stage for each cancer type; estimates were adjusted by isotonic regression weighted by the number of observations available per stage.'),
                  dataTableOutput("sens_table"),
                  p(HTML('Cancer: Cancer sensitivity group <br>
                   Stage: AJCC stage(I-IV) or NotStaged for no standard staging <br>
                   IR: Imputed incidence rate per 100K individuals from SEER <br>
                   Survival: 5 year cancer specific survival <br>
                   c: Cancers successfully detected <br>
                   n: Number of attempts <br>
                   sens: original sensitivity of multi-cancer early detection test <br>
                   iso_sens: Weighted isotonic regression estimate of sensitivity'))
                )
              )
      ),
      tabItem(tabName = "model_slip_parameters",
              fluidPage(
                box(
                  title = "Step 4: Set model slip parameters",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  fluidRow(
                    column(
                      6,
                      sliderInput("weibull_shape", label = h3("Weibull shape"), min = 0.001, max = 5, value = 1)
                    ),
                    column(
                      6,
                      sliderInput("screen_interval", label = h3("Screen interval (years)"), min = 0.5, max = 5, value = 1)
                    )
                  ),
                  plotOutput("weibull_dist_plot"),
                  p(HTML('Weibull shape: shape parameter of the Weibull distribution <br>
                   Screen interval: interval between screening events in years'))
                )
              )
      ),
      tabItem(tabName = "simulate_performance",
              fluidPage(
                box(
                  title = "Step 5: Simulate Performance",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  style = "overflow-x: auto;",
                  downloadButton("downloadTable4", "Download Table 4"),
                  dataTableOutput("simulated_data"),
                  p(HTML('Cancer: Cancer sensitivity group <br>
             clinical: Number of stage of original clinical presentation (1-4), 0 for NotStaged <br>
             prequel: Number of stage at interception (1-4), 0 for NotStaged <br>
             found_clinical: 1 = intercepted, 2 = clinical presentation <br>
             caught: Incidence intercepted at this prequel from clinical with method found_clinical <br>
             s_survival: shifted survival after interception <br>
             c_survival: original cancer survival <br>
             original_survivors: incidence surviving 5 years based on original stage <br>
             shifted_survivors: incidence surviving 5 years based on shifted stage <br>
             original_deaths: incidence dying after 5 years based on original stage <br>
             shifted_deaths: incidence dying after 5 years based on shifted stage <br>
             dw_scenario: dwell time scenario, "MIS": maximum interception scenario,"VSlow": very slow,"AggFast": aggressively fast <br>
             scan: Type of screening year: incident/prevalent/no screening <br>
             mode_found: cfdna or soc (matches found_clinical) cfdna: incidence found by cfdna screening; soc: incidence found by usual care<br>
             aggressive: dwell time scenario in words'))
                )
              )
      ),
      tabItem(tabName = "flow_diagrams",
              fluidPage(
                box(
                  title = "Step 6: Reconstruct Flow Diagrams",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  plotOutput("flow_diagram", height = "auto")
                )
              )
      ),
      tabItem(tabName = "utility_sankey_plot",
              fluidPage(
                box(
                  title = "Step 7: Utility Sankey Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  selectInput("sankey_plot_scenario", label = h3("Select scenario"), 
                              choices = list("MIS" = 1, "Very Slow" = 2, "Slow" = 3, 'Fast' = 4, 'Aggressively fast' = 5), 
                              selected = 1),
                  plotOutput("utility_sankey_plot", height = "auto")
                )
              )
      ),
      tabItem(tabName = "sensitivity_analysis",
              fluidPage(
                box(
                  title = "Step 8: Sensitivity Analysis by Cancer",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  downloadButton("download_stage_shift_table", "Download stage shift table"),
                  dataTableOutput("stage_shifting_table")
                )
              )
      )
    )
  )
)
