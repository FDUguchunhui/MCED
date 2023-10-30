library(shiny)
library(DT)
library(tidyverse)

source("../../scripts/current_date_code.R")
source("../../R/slip_rate_from_dwell.R")
source("../../R/run_intercept_model.R")

origin_host_dir = '../../.'
#read in standard performance numbers
base_detect_code<-"CCGA2"
#load and clean manuscript sensitivity
iso_sens_joined<-read_tsv(sprintf("%s/generated_data/%s_iso_seer_manuscript.tsv",origin_host_dir, input_date_code))
#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")



ui <- fluidPage(
  
  # Application title
  h1("Interception model"),
  p("This page will guide you through the process of simulated intercepted
    cancer cases and deaths at each stage and cancer type."),
  
  h1('Step1: Load dwell standard_model configuration csv file'),
  titlePanel("Load dwell standard_model configuration csv file"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose a CSV file"),
      br(),
      actionButton("loadFile", "Load File")
    ),
    
    mainPanel(
      DTOutput("table")
    )
  ),
  
  h1('Step 2: Isotonic regression adjusted Senstivity table'),
  dataTableOutput("sens_table"),
  
  h2("Dictionary of variables"),
  p(HTML('Cancer:	Cancer sensitivity group <br>
    Stage:	AJCC stage(I-IV) or NotStaged for no standard staging <br>
    IR:	Imputed incidence rate per 100K individuals from SEER <br>
    Survival:	5 year cancer specific survival <br>
    c:	Cancers successfully detected <br>
    n: Number of attempts <br>
    sens: original sensitivity of multi-cancer early dectection test <br>
    iso_sens:	Weighted isotonic regression estimate of sensitivity')),
  
  
  h1('Step 3: Set dwell time parameters'),
  
    
  fluidRow(
           # Copy the line below to make a slider bar 
           column(6, 
                  sliderInput("weibull_shape", label = h3("Weibull shape"), min = 0.001, 
                       max = 5, value = 1)
           ),
           
           column(6, 
                  sliderInput("screen_interval", label = h3("Screen interval"), min = 1, 
                              max = 10, value = 1))
  ),
  
  h1('Step 4: Simulate performance'),
  dataTableOutput("simulated_data"),
  h2("Dictionary of variables"),
  p(HTML('Cancer:	Cancer sensitivity group <br>
    clinical:	Number of stage of original clinical presentation (1-4), 0 for NotStaged <br>
    prequel:	Number of stage at interception (1-4), 0 for NotStaged <br>
    found_clinical:	1 = intercepted, 2 = clinical presentation <br>
    caught:	Incidence intercepted at this prequel from clinical with method found_clinical <br>
    s_survival: shifted survival after interception <br>
    c_survival: original cancer survival <br>
    original_survivors: incidence surviving 5 years based on original stage <br>
    shifted_survivors: incidence surviving 5 years based on shifted stage <br>
    original_deaths: incidence dying after 5 years based on original stage <br>
    shifted_deaths: incidence dying after 5 years based on shifted stage <br>
    dw_scenario: dwell time scenario <br>
    scan: Type of screening year: incident/prevalent/no screening <br>
    mode_found: cfdna or soc (matches found_clinical) <br>
    aggressive:	dwell time scenario in words')),
)

server <- function(input, output) {
  dwell_standard_model <- reactive({
    req(input$loadFile)
    infile <- input$file
    if (is.null(infile))
      return(NULL)
    
    read.csv(infile$datapath)
  })
  
  output$table <- renderDT({
    datatable(dwell_standard_model())
  })
  
 
  output$sens_table <- renderDataTable({
    incidence_sens_source
  })
  
  
  # add widgets control screen interval and weibull_shape or failure time distribution
  dwell_slip_rate <- reactive({
    exact_slip_rate_from_dwell(dwell_standard_model(),screen_interval=input$screen_interval,weibull_shape=input$weibull_shape)
  })
  
  
  #generate prevalent slip rate by clever use of very large interval and multiplying expectation
  long_interval<-100
  dwell_prevalent_rate <- reactive({
    exact_slip_rate_from_dwell(dwell_standard_model(),screen_interval=long_interval,weibull_shape=input$weibull_shape)
  })
  
  #no screening is happening, therefore nothing is ever intercepted
  dwell_no_rate <- reactive({
    dwell_slip_rate() %>% 
    mutate(slip=1.0,
           slip_clinical=1.0)
  })
  
  
  # generate simulated performance
  scenerios <- reactive({
    
    dwell_slip_rate <- dwell_slip_rate()
    dwell_no_rate <- dwell_no_rate()
    dwell_prevalent_rate <- dwell_prevalent_rate()
  #accumulate 4 dwell scenarios
  # prevalent and incident results
  # plus perfect screening
  # plus no screening
  my_list<-vector("list",4*2+2)
  k<-1
  for (dw_scenario in 1:4){
    print(k)
    local_performance<-run_intercept_model(incidence_sens_source,
                                           dwell_slip_rate %>% 
                                             filter(scenario==dw_scenario))
    local_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
    local_performance<-bind_rows(local_performance,local_excluded)
    
    incident_performance<-local_performance
    local_performance<-local_performance %>%
      mutate(screen_interval=input$screen_interval,
             dw_scenario=dw_scenario,
             scan="incident")
    
    my_list[[k]]<-local_performance
    k<-k+1
    
    #generate prevalent round starting off screening
    #only going to use caught by cfdna and combine with incident
    #because expected rates = average over all years, we can reverse identity
    #to obtain first-year screen by multiplying
    #rather than doing a special integral for prevalent rounds
    prevalent_performance<-run_intercept_model(incidence_sens_source,
                                               dwell_prevalent_rate %>% 
                                                 filter(scenario==dw_scenario))
    
    prevalent_performance<-prevalent_performance %>% 
      filter(found_clinical==1) %>%
      mutate(caught=caught*long_interval,
             original_survivors=original_survivors*long_interval,
             shifted_survivors=shifted_survivors*long_interval,
             original_deaths=original_deaths*long_interval,
             shifted_deaths=shifted_deaths*long_interval) %>%
      bind_rows(incident_performance %>% 
                  filter(found_clinical==2))
    
    prevalent_performance<-prevalent_performance %>%
      mutate(screen_interval=input$screen_interval,
             dw_scenario=dw_scenario,
             scan="prevalent")
    my_list[[k]]<-prevalent_performance
    k<-k+1
  }
  
  #this is the MIS scenario where schedule sensitivity is perfect so slip rates are 0
  optimal_performance<-run_intercept_model(incidence_sens_source,
                                           dwell_slip_rate %>% 
                                             filter(scenario==0))
  optimal_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
  optimal_performance<-bind_rows(optimal_performance,optimal_excluded)
  optimal_performance<-optimal_performance %>% mutate(opt="0",dw_scenario=0,scan="incident")
  
  #no screening so nothing found by cfdna operations
  no_screening_performance<-run_intercept_model(incidence_sens_source,
                                                dwell_no_rate %>%
                                                  filter(scenario==0))
  no_screening_performance<-bind_rows(no_screening_performance, optimal_excluded) %>%
    mutate(opt="NO",dw_scenario=0,scan="no")
  
  all_options_df<-bind_rows(my_list,.id="opt") %>%
    bind_rows(optimal_performance) %>%
    bind_rows(no_screening_performance)
  
  
  #Now we have the full, detailed data frame
  #add some helper text fields to clarify states represented by each line of the file 
  text_options_df<-all_options_df %>%
    select(opt,Cancer,clinical,prequel,found_clinical,
           caught,s_survival,c_survival, 
           original_survivors,shifted_survivors,
           original_deaths,shifted_deaths,
           screen_interval,dw_scenario,scan) %>%
    mutate(mode_found=c("cfdna","soc")[found_clinical],
           aggressive=c("MIS","VSlow","Slow","Fast","AggFast")[dw_scenario+1])
  
  return(text_options_df)
  })
  
  
  output$simulated_data <- renderDataTable({
    scenerios()
  })

}

shinyApp(ui, server)
