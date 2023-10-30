library(shiny)
library(DT)
library(tidyverse)

source("../../scripts/current_date_code.R")
source("../../R/slip_rate_from_dwell.R")
source("../../R/run_intercept_model.R")
source("../../R/reconstruct_flow_in_detail.R")
source("../../R/plot_flow_diagrams.R")


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
  p('Assumption: ensitivity estimates are nondecreasing by increasing stage for each cancer type;
  estimates were adjusted by isotonic regression weighted by the number of observations available per stage.'),
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
                  sliderInput("screen_interval", label = h3("Screen interval (years)"), min = 0.5, 
                              max = 10, value = 1))
  ),
  
  # plot the  Weibull distribution and plot failure probability based on input Weibull shape and screen interval
  plotOutput("weibull_dist_plot"),
  
  h1('Step 4: Simulate performance'),
  downloadButton("downloadTable4", "Download Table 4"),
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
  
  h1('Step 5: Reconstruct flow diagrams'),
  plotOutput("flow_diagram"),
)

server <- function(input, output) {
  
  dwell_standard_model <- reactive({
    # Default file path
    default_filepath <- "/Users/cgu3/Downloads/dwell_model_all_df (1).csv"
    
    # Check if the user has uploaded a new file
    
    infile <- input$file
    # If user has uploaded a file, use that; otherwise, use the default file
    if (is.null(infile$datapath)) {
      return(read.csv(default_filepath))
    } else {
      return(read.csv(infile$datapath))
    }
  })
  
  
  output$table <- renderDT({
    datatable(dwell_standard_model())
  })
  
 
  output$sens_table <- renderDataTable({
    incidence_sens_source
  })
  
  
  # plot the  Weibull distribution and plot failure probability based on input Weibull shape and screen interval
  output$weibull_dist_plot <- renderPlot({
    weibull_shape<-input$weibull_shape
    screen_interval<-input$screen_interval
    x<-seq(0,5,0.1)
    y<-dweibull(x,shape=weibull_shape,scale=screen_interval)
    plot(x,y,type="l",xlab="years",ylab="failure probability")
    abline(v=screen_interval,col="red")
    abline(h=0.5,col="red")
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
  
  
  
  # Step 4: Simulate performance
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
  
  # download the simulated data
  output$downloadTable4 <- downloadHandler(
    filename = function() {
      "text_data_set.csv"
    },
    content = function(file) {
      write.csv(scenerios(), file, row.names = FALSE)
    }
  )
  
  
  
  
  # plot flow diagram Figure 1 in paper
  plot_flow_diagram <- reactive({
    
    global_figure_scale<-7.5
    global_box_scale<-0.045
    global_text_scale<-1.35
    color_caught=c("intercept"="plum","clinical"="grey80")
    
    slip_rate_df <- dwell_slip_rate()
    
    diagram_setup<-set_up_diagram_clean(color_caught=color_caught)
    par(mfrow=c(2,2))
    # par(mar = c(0, 0, 1, 1))  # Adjust the margin (bottom, left, top, right)
    # par(oma = c(0, 0, 2, 0))  # Adjust the outer margin (bottom, left, top, right)
    
    #1A
    my_locus<-blank_intercept_para(diagram_setup,"State-Transition Graph",local_box_size=global_box_scale)
    
    #simply fill in the boxes with the identity
    text(diagram_setup$map$x,diagram_setup$map$y,diagram_setup$map$name,cex=global_text_scale)
    
    #1B
    dw_scenario<-0  #start with pure scenario
    
    local_slip_rate_df <- slip_rate_df %>% 
      filter(scenario==dw_scenario)
    
    my_title<-"No Interception"
    no_intercept_flow <- intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                           intercept_start_at_stage=0)
    #start plotting
    my_locus <- blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
    plot_object_flow_tuned(no_intercept_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
    
    
    ##1C
    dw_scenario<-0  #start with pure scenario
    
    local_slip_rate_df<-slip_rate_df %>% 
      filter(scenario==dw_scenario)
    
    my_title<-"Interception Model"
    mis_flow<-intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                  intercept_start_at_stage=4)
    my_locus<-blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
    plot_object_flow_tuned(mis_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
    
    ## 1D
    flow_up_to_stage<-4
    my_title<-"Interception Model: Fast"
    #now use a finite slip rate scenario
    dw_scenario<-3  #start with pure scenario
    local_slip_rate_df<-slip_rate_df %>% 
      filter(scenario==dw_scenario)
    
    fast_flow<-intercept_with_flow(incidence_sens_source,local_slip_rate_df,
                                   intercept_start_at_stage=4)
    #start plotting
    #start plotting
    my_locus<-blank_intercept_para(diagram_setup,my_title,local_box_size=global_box_scale)
    final_plot <- plot_object_flow_tuned(fast_flow,diagram_setup,my_locus,flow_up_to_stage=4,cex_scale=global_text_scale)
    return(final_plot)
  })
  
  output$flow_diagram <- renderPlot({
    plot_flow_diagram()
  }, width = 1000, height = 800)

}

shinyApp(ui, server)
