library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(ggalluvial)
library(ggplot2)
library(patchwork)

# source R files with local argument set to TRUE otherwise function from sourced file cannot access shiny app environment
source("../../scripts/current_date_code.R", local = TRUE)
source("../../R/slip_rate_from_dwell.R", local = TRUE)
source("../../R/run_intercept_model.R", local = TRUE)
source("../../R/reconstruct_flow_in_detail.R", local = TRUE)
source("../../R/plot_flow_diagrams.R", local = TRUE)
source("../../R/utility_sankey_plot.R", local = TRUE)

# global variables -----------------------------------------------------------
origin_host_dir = '../../.'
# read in standard performance numbers
base_detect_code <- "CCGA2"
# load and clean manuscript sensitivity
iso_sens_joined <- read_tsv(sprintf("%s/generated_data/%s_iso_seer_manuscript.tsv", origin_host_dir, input_date_code))
# remove everything not staged for initial analysis
incidence_sens_source <- iso_sens_joined %>% filter(Stage != "NotStaged") %>% mutate(iso_sens = sens)
# keep not staged for adding back
incidence_excluded_source <- iso_sens_joined %>% filter(Stage == "NotStaged")

# ui.R ----------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Interception Model"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidPage(
      box(
        title = "Introduction",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        p("This page will guide you through the process of simulated intercepted cancer cases and deaths at each stage and cancer type.")
      ),
      box(
        title = "Step 1: Load Dwell Standard Model Configuration CSV File",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        sidebarLayout(
          sidebarPanel(
            width = 3,
            fileInput("file", "Choose a CSV file"),
            br(),
            actionButton("loadFile", "Load File")
          ),
          mainPanel(DTOutput("table"))
        ),
        p(HTML('The dwell time for each cancer is defined by which dwell group it belong to and the progress scenario
               this table is to show the final configuration of the dwell time for each cancer type. <br>')),
        p(HTML('Cancer: Cancer sensitivity group <br>
               dwell_group: the pre-defined dwell group, each dwell group has a fixed dwell time for stage 1-4 <br>
               scenario: 1."VSlow": very slow, 2. slow 3. Fast 4. "AggFast": aggressively fast  <br>
               stage: the character stage name <br>
               number_stage: the numerical representative of the stage I: 1, II: 2, III: 3, IV: 4 <br>'))
      ),
      box(
        title = "Step 2: Isotonic Regression Adjusted Sensitivity Table",
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
      ),
      box(
        title = "Step 3: Set Dwell Time Parameters",
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
      ),
      box(
        title = "Step 4: Simulate Performance",
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
      ),
      box(
        title = "Step 5: Reconstruct Flow Diagrams",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        plotOutput("flow_diagram", height = "auto")
      ),
      box(
        title = "Step 6: Utility Sankey Plot",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        selectInput("sankey_plot_scenario", label = h3("Select scenario"), 
                    choices = list("MIS" = 1, "Very Slow" = 2, "Slow" = 3, 'Fast' = 4, 'Aggressively fast' = 5), 
                    selected = 1),
        plotOutput("utility_sankey_plot", height = "auto")
      ),
      box(
        title = "Step 7: Sensitivity Analysis by Cancer",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        downloadButton("download_stage_shift_table", "Download stage shift table"),
        dataTableOutput("stage_shifting_table")
      )
    )
  )
)

# server.R ------------------------------------------------------------------
server <- function(input, output) {

  
  # dwell standard model table reactive ------------------------------------------------
  dwell_standard_model <- reactive({
    # Default file path
    default_filepath <-
      "/Users/cgu3/Downloads/dwell_model_all_df (1).csv"
    
    # Check if the user has uploaded a new file
    
    infile <- input$file
    # If user has uploaded a file, use that; otherwise, use the default file
    if (is.null(infile$datapath)) {
      return(read.csv(default_filepath))
    } else {
      return(read.csv(infile$datapath))
    }
  })
  
  
  
  # Step1: output the dwell standard model table ------------------------------------------------
  output$table <- renderDT({
    datatable(dwell_standard_model())
  })
  
  # Step2: output sensitivity table ------------------------------------------------
  output$sens_table <- renderDataTable({
    incidence_sens_source
  })
  
  
  # plot the  Weibull distribution and plot failure probability based on input Weibull shape and screen interval ------------------------------------------------
  output$weibull_dist_plot <- renderPlot({
    par(mfrow=c(1,2))
    weibull_shape <- input$weibull_shape
    screen_interval <- input$screen_interval
    dwell_scale <- 4/gamma(1+1/weibull_shape)
    x <- seq(0, 5, 0.01)
    
    y_density <- dweibull(x, shape = weibull_shape, scale = dwell_scale)
    cum_prob <- pweibull(x, shape = weibull_shape, scale = dwell_scale)
    # Create data frames for ggplot
    density_data <- data.frame(x, y_density)
    cum_prob_data <- data.frame(x, cum_prob)
    
    base_theme <- theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)
      )
    
    # First plot: Weibull density
    p1 <- ggplot(density_data, aes(x = x, y = y_density)) +
      geom_line() +
      xlab("Years") +
      ylab("Slip density") +
      ggtitle("Slip density", subtitle = 'Plot with Dwell=4 for example') +
      base_theme
    
    # Second plot: Weibull cumulative probability
    p2 <- ggplot(cum_prob_data, aes(x = x, y = cum_prob)) +
      geom_line() +
      geom_vline(xintercept = screen_interval, color = "red", linetype = "dashed") +
      xlab("Years") +
      ylab("Cumulative Slip Probability") +
      ggtitle("Cumulative Slip Probability") +
      base_theme
    
    # Arrange the plots side by side using patchwork
    p1 + p2
  })
  
  
  # Step 3: Set dwell time parameters ------------------------------------------------
  ## dwell slip rate reactive ------------------------------------------------
  dwell_slip_rate <- reactive({
    exact_slip_rate_from_dwell(
      dwell_standard_model(),
      screen_interval = input$screen_interval,
      weibull_shape = input$weibull_shape
    )
  })
  
  ## dwell prevalent rate reactive ------------------------------------------------
  #generate prevalent slip rate by clever use of very large interval and multiplying expectation
  long_interval <- 100
  dwell_prevalent_rate <- reactive({
    exact_slip_rate_from_dwell(
      dwell_standard_model(),
      screen_interval = long_interval,
      weibull_shape = input$weibull_shape
    )
  })
  
  ## dwell no rate reactive ------------------------------------------------
  #no screening is happening, therefore nothing is ever intercepted
  dwell_no_rate <- reactive({
    dwell_slip_rate() %>%
      mutate(slip = 1.0,
             slip_clinical = 1.0)
  })
  
  
  
  # Step 4: Simulate performance ------------------------------------------------
  rich_survival <- reactive({
    dwell_slip_rate <- dwell_slip_rate()
    dwell_no_rate <- dwell_no_rate()
    dwell_prevalent_rate <- dwell_prevalent_rate()
    #accumulate 4 dwell scenarios
    # prevalent and incident results
    # plus perfect screening
    # plus no screening
    my_list <- vector("list", 4 * 2 + 2)
    k <- 1
    for (dw_scenario in 1:4) {
      print(k)
      
      #browser()
      
      local_performance <- run_intercept_model(incidence_sens_source,
                                               dwell_slip_rate %>%
                                                 filter(scenario == dw_scenario))
      local_excluded <-
        run_excluded_model(incidence_excluded_source) #does not depend on scenario
      local_performance <- bind_rows(local_performance, local_excluded)
      
      incident_performance <- local_performance
      local_performance <- local_performance %>%
        mutate(
          screen_interval = input$screen_interval,
          dw_scenario = dw_scenario,
          scan = "incident"
        )
      
      my_list[[k]] <- local_performance
      k <- k + 1
      
      #generate prevalent round starting off screening
      #only going to use caught by cfdna and combine with incident
      #because expected rates = average over all years, we can reverse identity
      #to obtain first-year screen by multiplying
      #rather than doing a special integral for prevalent rounds
      prevalent_performance <-
        run_intercept_model(
          incidence_sens_source,
          dwell_prevalent_rate %>%
            filter(scenario == dw_scenario)
        )
      
      prevalent_performance <- prevalent_performance %>%
        filter(found_clinical == 1) %>%
        mutate(
          caught = caught * long_interval,
          original_survivors = original_survivors * long_interval,
          shifted_survivors = shifted_survivors * long_interval,
          original_deaths = original_deaths * long_interval,
          shifted_deaths = shifted_deaths * long_interval
        ) %>%
        bind_rows(incident_performance %>%
                    filter(found_clinical == 2))
      
      prevalent_performance <- prevalent_performance %>%
        mutate(
          screen_interval = input$screen_interval,
          dw_scenario = dw_scenario,
          scan = "prevalent"
        )
      my_list[[k]] <- prevalent_performance
      k <- k + 1
    }
    
    #browser()
    
    #this is the MIS scenario where schedule sensitivity is perfect so slip rates are 0
    optimal_performance <- run_intercept_model(incidence_sens_source,
                                               dwell_slip_rate %>%
                                                 filter(scenario == 0))
    optimal_excluded <-
      run_excluded_model(incidence_excluded_source) #does not depend on scenario
    optimal_performance <-
      bind_rows(optimal_performance, optimal_excluded)
    optimal_performance <-
      optimal_performance %>% mutate(opt = "0",
                                     dw_scenario = 0,
                                     scan = "incident")
    
    #no screening so nothing found by cfdna operations
    no_screening_performance <-
      run_intercept_model(incidence_sens_source,
                          dwell_no_rate %>%
                            filter(scenario == 0))
    no_screening_performance <-
      bind_rows(no_screening_performance, optimal_excluded) %>%
      mutate(opt = "NO",
             dw_scenario = 0,
             scan = "no")
    
    all_options_df <- bind_rows(my_list, .id = "opt") %>%
      bind_rows(optimal_performance) %>%
      bind_rows(no_screening_performance)
    
    
    #Now we have the full, detailed data frame
    #add some helper text fields to clarify states represented by each line of the file
    text_options_df <- all_options_df %>%
      select(
        opt,
        Cancer,
        clinical,
        prequel,
        found_clinical,
        caught,
        s_survival,
        c_survival,
        original_survivors,
        shifted_survivors,
        original_deaths,
        shifted_deaths,
        screen_interval,
        dw_scenario,
        scan
      ) %>%
      mutate(
        mode_found = c("cfdna", "soc")[found_clinical],
        aggressive = c("MIS", "VSlow", "Slow", "Fast", "AggFast")[dw_scenario +
                                                                    1]
      )
    
    return(text_options_df)
  })
  
  
  output$simulated_data <- renderDataTable({
    rich_survival()
  }, options = list(pageLength = 10, scrollX = TRUE, scrollY = TRUE))
  
  # download the simulated data
  output$downloadTable4 <- downloadHandler(
    filename = function() {
      "text_data_set.csv"
    },
    content = function(file) {
      write.csv(rich_survival(), file, row.names = FALSE)
    }
  )
  
  
  

  # Step 5: plot flow diagram Figure 1 in paper -----------------------------------------
  
  plot_flow_diagram <- reactive({
    
    global_figure_scale <- 7.5
    global_box_scale <- 0.045
    global_text_scale <- 1.35
    color_caught <- c("intercept" = "plum", "clinical" = "grey80")
    
    slip_rate_df <- dwell_slip_rate()
    
    diagram_setup <- set_up_diagram_clean(color_caught = color_caught)
    par(mfrow = c(2, 2))
    
    # 1A
    my_locus <- blank_intercept_para(diagram_setup, "State-Transition Graph", local_box_size = global_box_scale)
    
    # Simply fill in the boxes with the identity
    text(diagram_setup$map$x, diagram_setup$map$y, diagram_setup$map$name, cex = global_text_scale)
    
    # 1B
    dw_scenario <- 0  # Start with pure scenario
    
    local_slip_rate_df <- slip_rate_df %>% filter(scenario == dw_scenario)
    
    my_title <- "No Interception"
    no_intercept_flow <- intercept_with_flow(incidence_sens_source, local_slip_rate_df, intercept_start_at_stage = 0)
    
    # Start plotting
    my_locus <- blank_intercept_para(diagram_setup, my_title, local_box_size = global_box_scale)
    plot_object_flow_tuned(no_intercept_flow, diagram_setup, my_locus, flow_up_to_stage = 4, cex_scale = global_text_scale)
    
    # 1C
    dw_scenario <- 0  # Start with pure scenario
    
    local_slip_rate_df <- slip_rate_df %>% filter(scenario == dw_scenario)
    
    my_title <- "Interception Model: MIS"
    mis_flow <- intercept_with_flow(incidence_sens_source, local_slip_rate_df, intercept_start_at_stage = 4)
    
    my_locus <- blank_intercept_para(diagram_setup, my_title, local_box_size = global_box_scale)
    plot_object_flow_tuned(mis_flow, diagram_setup, my_locus, flow_up_to_stage = 4, cex_scale = global_text_scale)
    
    # 1D
    flow_up_to_stage <- 4
    my_title <- "Interception Model: Fast"
    
    # Now use a finite slip rate scenario
    dw_scenario <- 3  # Start with pure scenario
    local_slip_rate_df <- slip_rate_df %>% filter(scenario == dw_scenario)
    
    fast_flow <- intercept_with_flow(incidence_sens_source, local_slip_rate_df, intercept_start_at_stage = 4)
    
    # Start plotting
    my_locus <- blank_intercept_para(diagram_setup, my_title, local_box_size = global_box_scale)
    final_plot <- plot_object_flow_tuned(fast_flow, diagram_setup, my_locus, flow_up_to_stage = 4, cex_scale = global_text_scale)
    
    return(final_plot)
  })
  
  output$flow_diagram <- renderPlot({
    plot_flow_diagram()
  }, width = 1000, height = 800)
  
  
  
  # Step 6: plot flow diagram Figure 3 in paper -----------------------------------------
  
  a_intercept <-reactive({
    text_levels <- c("MIS", "VSlow", "Slow", "Fast", "AggFast")
    scenario_code <-text_levels[as.integer(input$sankey_plot_scenario)]
    rich_survival() %>% filter(scan == "incident", aggressive == scenario_code)
  })

    
  
  plot_sankey_plot <- reactive({
    cStage <- c("NS", "I", "II", "III", "IV")
    text_levels <- c("MIS", "VSlow", "Slow", "Fast", "AggFast")
  
    scenario_code <-text_levels[as.integer(input$sankey_plot_scenario)]
    
    # a_intercept contains intercepted cases at each stage and cancer of all usual care diagnosed cancers during incidence round
    a_intercept <- a_intercept()
    
    #stage-shift
    intercept_shifted <- a_intercept %>%
      filter(mode_found == "cfdna") %>%
      group_by(Cancer, prequel) %>%
      summarize(caught = sum(caught)) %>%
      ungroup() %>%
      mutate(Stage = cStage[prequel + 1])
    
    intercept_original <- a_intercept %>%
      filter(mode_found == "cfdna") %>%
      group_by(Cancer, clinical) %>%
      summarize(caught = sum(caught)) %>%
      ungroup() %>%
      mutate(Stage = cStage[clinical + 1])
    
    #generate horizontal barplot
    
    stage_shift_in_intercepted_h_plot <- intercept_shifted %>%
      select(Cancer, Stage, cfdna = caught) %>%
      left_join(intercept_original %>%
                  select(Cancer, Stage, clinical = caught)) %>%
      gather(key = "case", value = "caught", cfdna, clinical) %>%
      mutate(case = case_when(case == "clinical" ~ "pre-intercept",
                              TRUE ~ "intercepted")) %>%
      group_by(case, Stage) %>%
      summarize(caught = sum(caught)) %>%
      ungroup() %>%
      mutate(Stage = factor(Stage, levels = rev(cStage))) %>%
      ggplot(aes(x = Stage, fill = case, y = caught)) +
      geom_col(width = 0.9,
               position = "dodge",
               color = "black") +
      geom_text(
        aes(label = round(caught, 0)),
        position = position_dodge(width = 0.9),
        hjust = -0.1,
        color = "black",
        size = 6
      ) +
      scale_fill_manual(values = c("plum", "grey80")) +
      scale_x_discrete(breaks = rev(cStage), drop = FALSE) +
      scale_y_continuous(position = "right") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        title = element_text(size = 18),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank()
      ) +
      labs(y = "Number of diagnoses per 100K") +
      coord_flip(ylim = c(0, 250)) +
      ggtitle(sprintf("%s: Intercepted with Stage Shift", scenario_code))
    
    stage_same_in_not_caught_h_plot <- a_intercept %>%
      filter(mode_found == "soc") %>%
      group_by(clinical) %>%
      summarize(caught = sum(caught)) %>%
      ungroup() %>%
      mutate(case = "clinical", Stage = cStage[clinical + 1]) %>%
      mutate(Stage = factor(Stage, levels = rev(cStage))) %>%
      mutate(case = "usual care") %>%
      ggplot(aes(x = Stage, fill = case, y = caught)) +
      geom_col(width = 0.9,
               position = "dodge",
               color = "black") +
      geom_text(
        aes(label = round(caught, 0)),
        position = position_dodge(width = 0.9),
        hjust = -0.1,
        color = "black",
        size = 6
      ) +
      scale_fill_manual(values = c("grey50")) +
      scale_x_discrete(breaks = rev(cStage), drop = FALSE) +
      scale_y_continuous(position = "right") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        title = element_text(size = 18),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank()
      ) +
      labs(y = "Number of diagnoses per 100K") +
      coord_flip(ylim = c(0, 325)) +
      ggtitle(sprintf("%s: Remaining Not Intercepted", scenario_code))
    
    #flip reverse
    joint_h_plot <-
      stage_shift_in_intercepted_h_plot + stage_same_in_not_caught_h_plot
    
    
    #reduce to stage found, stage at, mode found
    all_a_summary <- a_intercept %>%
      group_by(mode_found, clinical, prequel) %>%
      summarize(
        IR = sum(caught),
        Deaths = sum(shifted_deaths),
        Delta = sum(original_deaths),
        scan = scan[1],
        aggressive = aggressive[1]
      ) %>%
      mutate(Delta = Delta - Deaths,
             Survived = IR - Deaths - Delta) %>%
      ungroup()
    
    #generate some automatic scaling for option 0
    ee_plot_range <- all_a_summary %>%
      group_by(mode_found) %>%
      summarize(v = sum(IR)) %>%
      ungroup() %>%
      summarize(v = max(v)) %>%
      mutate(
        height_positive = round(v / 20 + 0.5, 0) * 20,
        height_ticks = pmin(100, round(pmax(
          height_positive / 4, 10
        ) / 10) * 10)
      )
    
    height_positive <- ee_plot_range$height_positive[1]
    height_ticks <- ee_plot_range$height_ticks[1]
    
    found_code <- "cfdna"
    all_b_summary <- all_a_summary %>% filter(mode_found == found_code)
    
    a_plot <-
      fancy_alluvial_plot_pct_more(
        all_b_summary,
        my_title = sprintf("%s: Intercepted With Mortality Shift", scenario_code),
        height_positive,
        height_ticks
      )
    
    
    found_code <- "soc"
    all_b_summary <- all_a_summary %>% filter(mode_found == found_code)
    
    b_plot <-
      fancy_alluvial_plot_pct_more(
        all_b_summary,
        my_title = sprintf("%s: Remaining With Mortality", scenario_code),
        height_positive,
        height_ticks
      )
    
    
    total_figure3_top <-
      (stage_shift_in_intercepted_h_plot |
         stage_same_in_not_caught_h_plot) + plot_layout(tag_level = 'new')
    total_figure3_bottom <-
      (a_plot | b_plot) + plot_layout(tag_level = 'new')
    total_figure3_plot <-
      total_figure3_top / total_figure3_bottom + plot_annotation(tag_levels =
                                                                   c('A', '1'))  &
      theme(plot.tag = element_text(size = 16))
    
    return(total_figure3_plot)
  })
  
  
  output$utility_sankey_plot <- renderPlot({
    plot_sankey_plot()
  },  width = 1600, height = 1000)
  
  
  
  # Step 7: display stage shifting table -----------------------------
  
  sens_table <- reactive({
    a_intercept() %>% 
      filter(Cancer %in% c('Breast', 'Colon/Rectum', 'Esophagus', 'Liver/Bile-duct', 'Lung',
                           'Ovary', 'Pancreas', 'Stomach')) %>% 
      group_by(Cancer, clinical, prequel) %>% summarise(caught=sum(caught)) %>%
      group_by(Cancer, clinical) %>%
      mutate(total_at_corresponding_prequel = sum(caught)) %>%
      mutate(prob_at_prequel = caught / total_at_corresponding_prequel) %>%
      mutate(prob_at_prequel_percentage = scales::percent(prob_at_prequel, accuracy = 0.01)) %>%  
      select(Cancer, clinical, prequel, prob_at_prequel_percentage) %>% arrange(Cancer, clinical, prequel) %>%
      rename(early_detection = prequel, prob = prob_at_prequel_percentage) %>% 
      pivot_wider(names_from = clinical, values_from = prob) %>%
      arrange(Cancer, early_detection) %>% 
      filter(!if_all(c("1", "2", "3", "4"), ~ .x == "")) %>% 
      rename(`Stage 1` = `1`, `Stage 2` = `2`, `Stage 3` = `3`, `Stage 4` = `4`)
    
  })
  
  output$stage_shifting_table <- renderDataTable({
    sens_table()
  })
  
  
  output$download_stage_shift_table <- downloadHandler(
    filename = function() {
      "stage_shift_table.csv"
    },
    content = function(file) {
      write.csv(sens_table(), file, row.names = FALSE)
    }
  )
  
  
  
}

shinyApp(ui, server)
