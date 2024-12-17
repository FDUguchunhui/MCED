xStage <- c("I", "II", "III", "IV")

# dictionary:
# prequel: stage at early-detection test
# clinical: clinical stage/original clinical presentation stage/clinically diagnosed stage/usual-care 
# dignosed stage

#okay: this multiplies a final destination IR: all the people who wind up at some final destination
#so everything scales within itself
#assume stage IV destination unless the cases are stopped earlier

#' @title Calculate effective detection rate with consideration of slip rate
#' @description This function takes in cumulative sensitivity `cumulative_sens` and dwell detect rates
#' `dwell_detect_rate` to calculate the detection rate at each stage and the number of cases missed.
#' @param cumulative_sens cumulative sensitivity
#' @param dwell_detect_rate dwell detect rate
#' @details
#' For example, when clinical stage = 4, prequel can be 1, 2, 3, 4
#' cumulative_sens
#' [1] 0.3333333 0.6000000 1.0000000 1.0000000
#' dwell_detect_rate
#' [1] 0.8847969 0.7869388 0.6321210 0.4323334
effective_sens <- function(cumulative_sens, dwell_detect_rate) {
  detect <- rep(0, length(cumulative_sens))
  miss <- 0
  i <- 1
  arrive <- cumulative_sens[i]
  live <- arrive + miss
  detect[1] <- cumulative_sens[i] * dwell_detect_rate[i]
  miss <- cumulative_sens[i] * (1 - dwell_detect_rate[i])
  if (length(cumulative_sens) > 1) {
    for (i in 2:length(cumulative_sens))
    {
      #newly detectable cases
      arrive <- cumulative_sens[i] - cumulative_sens[i - 1]
      live <-
        arrive + miss #currently detectable is newly detectable + missed at earlier stages
      detect[i] <-
        live * dwell_detect_rate[i] #would detect all of them, but miss some of them due to timing
      miss <- live * (1 - dwell_detect_rate[i])
    }
  }
  arrive <- 1 - cumulative_sens[i] #final miss
  live = arrive + miss
  #detect[1-4] is detected at each stage
  #clinical detect= miss[4]
  list(intercept = detect, clinical = live)
}


#' @title 
#' @description  It adds survival rates to the detection stages, using 5-year survival rates as a rough estimate of
#' "statistical cure" unaffected by 1-2 years of lead time (time between disease detectability and diagnosis).
#' It then computes absolute numbers of survivors and deaths before and after a hypothetical shift in detection stage
#' @param incidence_sens_source incidence_sens_source
#' @param incidence_intercepted incidence_intercepted
#' @details
#' #' #' incidence_sens_source
#' | Cancer   | Stage | IR    | Survival | c | n | sens  | iso_sens |
#' |----------|-------|-------|----------|---|---|-------|----------|
#' | Anus     | I     | 1.20  | 0.909    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | II    | 1.82  | 0.821    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | III   | 1.61  | 0.672    | 3 | 5 | 0.6   | 0.6      |
#' | Anus     | IV    | 0.503 | 0.235    | 0 | 0 | 1     | 1        |
#' | Bladder  | I     | 12.1  | 0.862    | 3 | 5 | 0.545 | 0.545    |
#' | Bladder  | II    | 5.50  | 0.618    | 2 | 3 | 0.545 | 0.545    |
#' | Bladder  | III   | 2.20  | 0.495    | 1 | 2 | 0.545 | 0.545    |
#' | Bladder  | IV    | 4.06  | 0.163    | 0 | 1 | 0.545 | 0.545    |
#' 
#' # incidence_intercepted:
# Cancer| Stage| clinical| prequel| detect  |delta_detect| modified_slip| sens_slip| found_clinical|  IR | caught
# ------|------| --------|--------|-------- |------------|--------------|----------|---------------|-----|-------
# Anus  | I    |  1      | 1      |  0.333  |    0.202   |   0.393      |  0.202   |    1          | 1.2 | 0.242
# Anus  | I    |  1      | 1      |  0.333  |    0.798   |   0.393      |  0.202   |    2          | 1.2 | 0.956
# Anus  | II   |  2      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 1.8 | 0.472
# Anus  | II   |  2      | 2      |  0.6    |    0.125   |   0.632      |  0.125   |    1          | 1.8 | 0.228
# Anus  | II   |  2      | 2      |  0.6    |    0.615   |   0.632      |  0.125   |    2          | 1.8 | 1.12
# Anus  | III  |  3      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 1.6 | 0.418
# Anus  | III  |  3      | 2      |  0.6    |    0.206   |   0.393      |  0.206   |    1          | 1.6 | 0.333
# Anus  | III  |  3      | 3      |  1      |    0.0723  |   0.865      |  0.0723  |    1          | 1.6 | 0.116
# Anus  | III  |  3      | 3      |  1      |    0.462   |   0.865      |  0.0723  |    2          | 1.6 | 0.744
# Anus  | IV   |  4      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 0.50| 0.130
add_survival_to_stage_shift <-
  function(incidence_sens_source,
           incidence_intercepted) {
    # browser()
    
    #add survival: original before shift, and shifted survival
    #using 5 year survival as "crude estimate of statistical cure"
    #unlikely to be affected by 1-2 year lead time significantly
    #extract survival by 'stage_at_detection'
   # just_survival:
   # |    |  Cancer  | prequel| Survival
   # | 1  |  Anus    |   1    |  0.909
   # | 2  |  Anus    |   2    |  0.821
   # | 3  |  Anus    |   3    |  0.672
   # | 4  |  Anus    |   4    |  0.235
   # | 5  |  Bladder |   1    |  0.862
   # | 6  |  Bladder |   2    |  0.618
   # | 7  |  Bladder |   3    |  0.495
   # | 8  |  Bladder |   4    |  0.163
   # | 9  |  Breast  |   1    |  0.982
   # | 10 |  Breast  |   2    |  0.929
    just_survival <- incidence_sens_source %>%
      mutate(prequel = match(Stage, xStage)) %>%
      select(Cancer, prequel, Survival) %>%
      filter(!is.na(prequel))
    
    #  intercept_survival:
    #   Cancer| Stage|    IR |number_stage| clinical| prequel| found_clinical| caught| s_survival| c_survival
    # 1  Anus |  I   |  1.20 |      1     |   1     |  1     |        1      |  0.242|      0.909|    0.909
    # 2  Anus |  I   |  1.20 |      1     |   1     |  1     |        2      |  0.956|      0.909|    0.909
    # 3  Anus |  II  |  1.82 |      2     |   2     |  1     |        1      |  0.472|      0.909|    0.821
    # 4  Anus |  II  |  1.82 |      2     |   2     |  2     |        1      |  0.228|      0.821|    0.821
    # 5  Anus |  II  |  1.82 |      2     |   2     |  2     |        2      |  1.12 |      0.821|    0.821
    # 6  Anus |  III |  1.61 |      3     |   3     |  1     |        1      |  0.418|      0.909|    0.672
    # 7  Anus |  III |  1.61 |      3     |   3     |  2     |        1      |  0.333|      0.821|    0.672
    # 8  Anus |  III |  1.61 |      3     |   3     |  3     |        1      |  0.116|      0.672|    0.672
    # 9  Anus |  III |  1.61 |      3     |   3     |  3     |        2      |  0.744|      0.672|    0.672
    # 10 Anus |  IV  | 0.503 |      4     |   4     |  1     |        1      |  0.130|      0.909|    0.235
    intercept_survival <- incidence_intercepted %>%
      left_join(
        just_survival %>%
          select(Cancer, prequel, s_survival = Survival),
        by = c("Cancer", "prequel")
      ) %>%
      left_join(
        just_survival %>%
          select(Cancer, clinical = prequel, c_survival = Survival),
        by = c("Cancer", "clinical")
      )
    
    #compute absolute numbers rather than local rates
    # five year survival rate and death rate (per 1 million)
    # intercept_survival:
    # Cancer | Stage |  IR  | caught | c_survival | s_survival | original_survivors | shifted_survivors | original_deaths| shifted_deaths
    # -------|-------|------|--------|------------|------------|--------------------|-------------------|----------------|--------------
    # Anus   |   I   | 1.20 |  0.242 |      0.909 |      0.909 |       0.220        |             0.220 |      0.0221    |    0.0221
    # Anus   |   I   | 1.20 |  0.956 |      0.909 |      0.909 |       0.869        |             0.869 |      0.0872    |    0.0872
    # Anus   |  II   | 1.82 |  0.472 |      0.821 |      0.909 |       0.388        |             0.429 |      0.0434    |    0.0434
    # Anus   |  II   | 1.82 |  0.228 |      0.821 |      0.821 |       0.187        |             0.187 |      0.0406    |    0.0406
    # Anus   |  II   | 1.82 |  1.12  |      0.821 |      0.821 |       0.918        |             0.918 |       0.200    |     0.200
    # Anus   | III   | 1.61 |  0.418 |      0.672 |      0.909 |       0.281        |             0.380 |       0.137    |    0.1381
    # Anus   | III   | 1.61 |  0.333 |      0.672 |      0.821 |       0.224        |             0.273 |       0.109    |    0.0594
    # Anus   | III   | 1.61 |  0.116 |      0.672 |      0.672 |       0.0783       |             0.0783|      0.0381    |    0.0381
    # Anus   | III   | 1.61 |  0.744 |      0.672 |      0.672 |       0.500        |             0.500 |       0.244    |     0.244
    # Anus   |  IV   | 0.503|  0.130 |      0.235 |      0.909 |       0.0307       |             0.119 |      0.0998    |    0.0119
    intercept_survival <- intercept_survival %>%
      mutate(
        original_survivors = c_survival * caught,
        shifted_survivors = s_survival * caught,
        original_deaths = (1 - c_survival) * caught,
        shifted_deaths = (1 - s_survival) * caught
      )
    intercept_survival
  }


#' @title
#' @description  This function calculates the effective detection rate by stage, taking into account the "slip rate" 
#' (how cases might slip through detection at each stage) and any adjustments for cases detected at the clinical stage.
#' @param incidence_sens_source incidence_sens_source
#' @param dwell_slip_df dwell_slip_df
#' @param active_slip_clinical active_slip_clinical
#' @details
#' 
#' #' incidence_sens_source
#' | Cancer   | Stage | IR    | Survival | c | n | sens  | iso_sens |
#' |----------|-------|-------|----------|---|---|-------|----------|
#' | Anus     | I     | 1.20  | 0.909    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | II    | 1.82  | 0.821    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | III   | 1.61  | 0.672    | 3 | 5 | 0.6   | 0.6      |
#' | Anus     | IV    | 0.503 | 0.235    | 0 | 0 | 1     | 1        |
#' | Bladder  | I     | 12.1  | 0.862    | 3 | 5 | 0.545 | 0.545    |
#' | Bladder  | II    | 5.50  | 0.618    | 2 | 3 | 0.545 | 0.545    |
#' | Bladder  | III   | 2.20  | 0.495    | 1 | 2 | 0.545 | 0.545    |
#' | Bladder  | IV    | 4.06  | 0.163    | 0 | 1 | 0.545 | 0.545    |
#' 
#' #' dwell_slip_df
#' | Cancer        |dwell_group | scenario | Stage | dwell | number_stage | slip  | slip_clinical | screen_interval |
#' |---------------|     ---    |   ---    |  -----|   --- |-      -      |-------|    -------  |       ---         |
#' | Anus          |      A     |    1     |   I   |    4  |       1      | 0.115 |     0.213   |        1          |
#' | Anus          |      A     |    1     |   II  |    2  |       2      | 0.213 |     0.368   |        1          |
#' | Anus          |      A     |    1     |   III |    1  |       3      | 0.368 |     0.568   |        1          |
#' | Anus          |      A     |    1     |   IV  |    1  |       4      | 0.368 |     0.568   |        1          |
#' | Colon/Rectum  |      A     |    1     |   I   |    4  |       1      | 0.115 |     0.213   |        1          |
#' | Colon/Rectum  |      A     |    1     |   II  |    2  |       2      | 0.213 |     0.368   |        1          |
#' | Colon/Rectum  |      A     |    1     |   III |    1  |       3      | 0.368 |     0.568   |        1          |
#' | Colon/Rectum  |      A     |    1     |   IV  |    1  |       4      | 0.368 |     0.568   |        1          |
#' 
#' @return a table of the cumulative detection rate at each stage
compute_effective_detection_with_slip <-
  function(incidence_sens_source,
           dwell_slip_df,
           active_slip_clinical) {
    
    
    
    #just detection rate
    # just_detection: is a table of the cumulative detection rate at each stage
    # | Cancer  | prequel | detect |
    #   |--------|---------|--------|
    #   | Anus   | 1       | 0.333  |
    #   | Anus   | 2       | 0.6    |
    #   | Anus   | 3       | 1      |
    #   | Anus   | 4       | 1      |
    #   | Bladder| 1       | 0.545  |
    #   | Bladder| 2       | 0.545  |
    #   | Bladder| 3       | 0.545  |
    #   | Bladder| 4       | 0.545  |
    just_detection <- incidence_sens_source %>%
      mutate(
        number_stage = match(Stage, xStage),
        prequel = number_stage,
        detect = iso_sens
      ) %>%
      select(Cancer, prequel, detect)
    
    #differences - marginal detection rate of remaining cases
    #given that cases already detectable at earlier stage were removed or treated separately
    # | Cancer  | prequel | detect | delta_detect |
    #   |--------|---------|--------|--------------|
    #   | Anus   | 1       | 0.333  | 0.333        |
    #   | Anus   | 2       | 0.6    | 0.267        |
    #   | Anus   | 3       | 1      | 0.4          |
    #   | Anus   | 4       | 1      | 0            |
    #   | Bladder| 1       | 0.545  | 0.545        |
    #   | Bladder| 2       | 0.545  | 0            |
    #   | Bladder| 3       | 0.545  | 0            |
    #   | Bladder| 4       | 0.545  | 0            |
    just_delta <- just_detection %>%
      group_by(Cancer) %>%
      arrange(prequel, .by_group = TRUE) %>%
      mutate(delta_detect = diff(c(0, detect))) %>%
      ungroup() %>%
      arrange(Cancer)
    
    # get sensitivity with slip
    # just_slip_delta_extra: is a table of the sensitivity with slip for at each prequel stage and possible clinical stage
    # example when `just_slip_delta_extra=TRUE`
    #  | Cancer | prequel | detect | delta_detect | slip  | slip_clinical | clinical | modified_slip | sens_slip |
    #   |--------|---------|--------|--------------|-------|---------------|----------|---------------|----------|
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 1        | 0.213         | 0.262    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 2        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 2        | 0.368         | 0.193    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 3        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 3        | 0.213         | 0.240    |
    #   | Anus   | 3       | 1      | 0.4          | 0.368 | 0.568         | 3        | 0.568         | 0.201    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 4        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 4        | 0.213         | 0.240    |
    #   | Anus   | 3       | 1      | 0.4          | 0.368 | 0.568         | 4        | 0.368         | 0.294    |
    #   | Anus   | 4       | 1      | 0            | 0.368 | 0.568         | 4        | 0.568         | 0.0740   |
    # This code is taking the just_delta dataframe and joining it with dwell_slip_df on the columns Cancer and prequel.
    # The dwell_slip_df dataframe is being transformed to only select the columns Cancer, number_stage, slip, 
    # and slip_clinical, and it's renaming the number_stage column to prequel
    # The code is adding a new column called unroll with a constant value of 4. Then use it to create 4 copy of each row
    # with a new column called `clinical` and only keep rows where clinical is less than prequel. `Clinical` is the 
    # clinical stage.
    just_slip_delta_extra <- just_delta %>%
      left_join(
        dwell_slip_df %>%
          #filter(scenario==dw_scenario) %>%
          select(Cancer, prequel = number_stage, slip, slip_clinical),
        by = c("Cancer", "prequel")
      ) %>%
      filter(!is.na(prequel)) %>%
      mutate(unroll = 4) %>%
      uncount(unroll, .id = "clinical") %>%
      filter(clinical >= prequel) %>%
      # a new column modified_slip is being added based on whether the prequel is less than clinical and whether the 
      # `active_slip_clinical` is true or not
      mutate(
        modified_slip = case_when(
          prequel < clinical ~ slip,
          prequel == clinical &
            active_slip_clinical ~ slip_clinical,
          prequel == clinical &
            !active_slip_clinical ~ slip,
          TRUE ~ 1.0
        )
      ) %>%
      arrange(Cancer, clinical, prequel) %>%
      group_by(Cancer, clinical) %>%
      #a new column sens_slip is created by applying the effective_sens function on the detect column and 
      #(1 - modified_slip). sens_slip is the effective sensitivity with slip
      mutate(sens_slip = effective_sens(detect, 1 - modified_slip)$intercept) %>%
      ungroup()
    
    just_slip_delta_extra
  }



#' @title
#' @description  This function calculates the effective detection rate by stage, taking into account the "slip rate"
#' 
#' @param incidence_sens_source incidence_sens_source
#' @param dwell_slip_df dwell_slip_df
#' 
#' @details
#' incidence_sens_source
#' | Cancer   | Stage | IR    | Survival | c | n | sens  | iso_sens |
#' |----------|-------|-------|----------|---|---|-------|----------|
#' | Anus     | I     | 1.20  | 0.909    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | II    | 1.82  | 0.821    | 1 | 3 | 0.333 | 0.333    |
#' | Anus     | III   | 1.61  | 0.672    | 3 | 5 | 0.6   | 0.6      |
#' | Anus     | IV    | 0.503 | 0.235    | 0 | 0 | 1     | 1        |
#' | Bladder  | I     | 12.1  | 0.862    | 3 | 5 | 0.545 | 0.545    |
#' | Bladder  | II    | 5.50  | 0.618    | 2 | 3 | 0.545 | 0.545    |
#' | Bladder  | III   | 2.20  | 0.495    | 1 | 2 | 0.545 | 0.545    |
#' | Bladder  | IV    | 4.06  | 0.163    | 0 | 1 | 0.545 | 0.545    |
  
#' 
#' dwell_slip_df
#' | Cancer        |dwell_group | scenario | Stage | dwell | number_stage | slip  | slip_clinical | screen_interval |
#' |---------------|     ---    |   ---    |  -----|   --- |-      -      |-------|    -------  |       ---         |
#' | Anus          |      A     |    1     |   I   |    4  |       1      | 0.115 |     0.213   |        1          |
#' | Anus          |      A     |    1     |   II  |    2  |       2      | 0.213 |     0.368   |        1          |
#' | Anus          |      A     |    1     |   III |    1  |       3      | 0.368 |     0.568   |        1          |
#' | Anus          |      A     |    1     |   IV  |    1  |       4      | 0.368 |     0.568   |        1          |
#' | Colon/Rectum  |      A     |    1     |   I   |    4  |       1      | 0.115 |     0.213   |        1          |
#' | Colon/Rectum  |      A     |    1     |   II  |    2  |       2      | 0.213 |     0.368   |        1          |
#' | Colon/Rectum  |      A     |    1     |   III |    1  |       3      | 0.368 |     0.568   |        1          |
#' | Colon/Rectum  |      A     |    1     |   IV  |    1  |       4      | 0.368 |     0.568   |        1          |
  
run_intercept_model <-
  function(incidence_sens_source,
           dwell_slip_df,
           active_slip_clinical = TRUE) {

    
    #set up all previous stages where cases could be intercepted given clinical detection
    # number_stage = match(Stage, xStage): For each value in the Stage column, match finds its position in the xStage 
    # vector or list. The result is stored in the new column number_stage.
    # unroll = number_stage: Duplicating the number_stage column and naming it unroll.
    # uncount(unroll, .id = "prequel"): The uncount function is used to expand or "blow up" the data frame by 
    # replicating each row based on the value in the unroll column, and create a new column "prequel" which lists all
    # possible previous stage before a given clinical stage in numeric values
    # the acquired incidence_set will be used to calculate intercepted cases by stage
    # incidence_set:
    # +--------+-------+-----+-------------+---------+--------+
    #   | Cancer | Stage | IR  | number_stage| clinical| prequel|
    #   +--------+-------+-----+-------------+---------+--------+
    #   | Anus   |   I   | 1.20|      1      |    1    |   1    |
    #   | Anus   |  II   | 1.82|      2      |    2    |   1    |
    #   | Anus   |  II   | 1.82|      2      |    2    |   2    |
    #   | Anus   | III   | 1.61|      3      |    3    |   1    |
    #   | Anus   | III   | 1.61|      3      |    3    |   2    |
    #   | Anus   | III   | 1.61|      3      |    3    |   3    |
    #   | Anus   |  IV   | 0.503|     4      |    4    |   1    |
    #   | Anus   |  IV   | 0.503|     4      |    4    |   2    |
    #   +--------+-------+-----+-------------+---------+--------+
    incidence_set <- incidence_sens_source %>%
      filter(Stage %in% xStage) %>%
      select(Cancer, Stage, IR) %>%
      mutate(
        number_stage = match(Stage, xStage),
        clinical = number_stage,
        unroll = number_stage
      ) %>%
      uncount(unroll, .id = "prequel")
    
    #compute effective detection by stage conditional on slip rate model
    # just_slip_delta_extra:
    #  | Cancer | prequel | detect | delta_detect | slip  | slip_clinical | clinical | modified_slip | sens_slip |
    #   |--------|---------|--------|--------------|-------|---------------|----------|---------------|----------|
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 1        | 0.213         | 0.262    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 2        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 2        | 0.368         | 0.193    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 3        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 3        | 0.213         | 0.240    |
    #   | Anus   | 3       | 1      | 0.4          | 0.368 | 0.568         | 3        | 0.568         | 0.201    |
    #   | Anus   | 1       | 0.333  | 0.333        | 0.115 | 0.213         | 4        | 0.115         | 0.295    |
    #   | Anus   | 2       | 0.6    | 0.267        | 0.213 | 0.368         | 4        | 0.213         | 0.240    |
    #   | Anus   | 3       | 1      | 0.4          | 0.368 | 0.568         | 4        | 0.368         | 0.294    |
    #   | Anus   | 4       | 1      | 0            | 0.368 | 0.568         | 4        | 0.568         | 0.0740   |
    just_slip_delta_extra <-
      compute_effective_detection_with_slip(incidence_sens_source,
                                            dwell_slip_df,
                                            active_slip_clinical)
    
  
    # browser()
    
    #updated: split effective slip rate in 2 for last stage as the "time spent" should be halved
    #this involves a more elaborate model
    #note that "lives saved" is not affected, because those individuals are not stage shifted
    #this assumes that 'stage 4' is just automatically halved anyway
    # incidence_intercepted:
    # Cancer| Stage| clinical| prequel| detect  |delta_detect| modified_slip| sens_slip| found_clinical|  IR | caught
    # ------|------| --------|--------|-------- |------------|--------------|----------|---------------|-----|-------
    # Anus  | I    |  1      | 1      |  0.333  |    0.202   |   0.393      |  0.202   |    1          | 1.2 | 0.242
    # Anus  | I    |  1      | 1      |  0.333  |    0.798   |   0.393      |  0.202   |    2          | 1.2 | 0.956
    # Anus  | II   |  2      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 1.8 | 0.472
    # Anus  | II   |  2      | 2      |  0.6    |    0.125   |   0.632      |  0.125   |    1          | 1.8 | 0.228
    # Anus  | II   |  2      | 2      |  0.6    |    0.615   |   0.632      |  0.125   |    2          | 1.8 | 1.12
    # Anus  | III  |  3      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 1.6 | 0.418
    # Anus  | III  |  3      | 2      |  0.6    |    0.206   |   0.393      |  0.206   |    1          | 1.6 | 0.333
    # Anus  | III  |  3      | 3      |  1      |    0.0723  |   0.865      |  0.0723  |    1          | 1.6 | 0.116
    # Anus  | III  |  3      | 3      |  1      |    0.462   |   0.865      |  0.0723  |    2          | 1.6 | 0.744
    # Anus  | IV   |  4      | 1      |  0.333  |    0.260   |   0.221      |  0.260   |    1          | 0.50| 0.130
                   
    incidence_intercepted <- incidence_set %>%
      left_join(just_slip_delta_extra, by = c("Cancer", "clinical", "prequel")) %>%
      # Duplicates rows when clinical stage == prequel into two, and the new found_clinical will take value 1 or 2.
      # 1 means the individual is detected by early-detection screen, 2 means the individual is found clinically
      # for no duplicated row, found_clinical will be 1 (found by early-detection screen)
      mutate(unroll = 1 + (number_stage == prequel)) %>%
      uncount(unroll, .id = "found_clinical") %>%
      group_by(Cancer, clinical) %>%
      arrange(Cancer, clinical, prequel, found_clinical) %>%
      # delta_detect is the new slip-adjusted probability of being detected
      mutate(
        c_slip = cumsum(sens_slip),
        # c_slip calculate cumulative slip-adjusted sensitivity for group (Cancer, clinical) for all prequel stages
        # since the last row of each group is the expanded row for clinical stage == prequel for remaining individuals
        # that diagnosed clinically at original clinical presentation stage, it has to be removed when calculating delta_detect in found_clinical == 2
        # so found_clinical == 2 ~ 1 - c_slip + sens_slip has a "plus sens_slip" part
        # delta_detect in this table represent percentage of cancers that are detected by early-detection screen at
        # each stage after adjusting for modified slip rate
        delta_detect = case_when(found_clinical == 2 ~ 1 - c_slip + sens_slip, #anyone not caught by new screening must be found clinically
                                 TRUE ~ sens_slip)
      ) %>%
      # caught is the number of individuals diagnosed in this prequel stage
      # when prequel == clinical individuals that haven't detected by early-detection screen will be found clinically
      mutate(caught = IR * delta_detect) %>%
      ungroup()
    
    intercept_survival <-
      add_survival_to_stage_shift(incidence_sens_source, incidence_intercepted)
    
    intercept_survival
  }


#' @title Fill in performance results for cancers without stage information
#' @description all cases are considered as clinical diagnosed (found_clinical = 2)
#' @param excluded_source data frame with individuals not staged as they are not modeled
#' @details
#' there is no survival and death improvement for using early-detection (s_survival = c_survival, s_death = c_death)
#' 
#' 
#' @return data frame with survival estimates
run_excluded_model <- function(excluded_source) {
  #fills out individuals not staged as they are not modeled
  excluded_survival <- excluded_source %>%
    mutate(
      number_stage = 0,
      clinical = 0,
      prequel = 0,
      detect = 0.0,
      delta_detect = 0.0,
      slip = 1.0,
      slip_clinical = 1.0,
      modified_slip = 1.0,
      sens_slip = 0.0,
      found_clinical = 2,
      c_slip = 1.0,
      caught = IR,
      s_survival = Survival,
      c_survival = Survival
    ) %>%
    mutate(
      original_survivors = c_survival * caught,
      shifted_survivors = s_survival * caught,
      original_deaths = (1 - c_survival) * caught,
      shifted_deaths = (1 - s_survival) * caught
    ) %>%
    select(
      Cancer,
      Stage,
      IR,
      number_stage,
      clinical,
      prequel,
      detect,
      delta_detect,
      slip,
      slip_clinical,
      modified_slip,
      sens_slip,
      found_clinical,
      c_slip,
      caught,
      s_survival,
      c_survival,
      original_survivors,
      shifted_survivors,
      original_deaths,
      shifted_deaths
    )
  
  excluded_survival
}
