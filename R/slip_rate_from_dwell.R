#slip rate computation

#exact: integrate over those missed during a screening interval
#see supplemental information for argument that this function is appropriate
#' @title Integrate slip rate
#' @description Integrate slip rate
#' @param screen_interval screening interval in years
#' @param weibull_shape weibull shape parameter
#' @param dwell dwell time in years
#' @details
#' 
integrate_slip_rate<-function(screen_interval, weibull_shape, dwell){
  # browser()
  # slip rate
  # yield of escape = integral cumululate distribution function F(t), 0<t<screen_interval
  # mean_of_weibull = weibull_scale*gamma(1+1/shape)
  # so weibull_scale = mean_of_weibull/gamma(1+1/shape) and mean_of_weibull = dwell, i.e. the average time of dwell/sojourn
  dwell_scale=dwell/gamma(1+1/weibull_shape)
  tiny_delta<-365
  #trapezoidal integration
  # https://www.onlinemathlearning.com/trapezoid-rule.html
  # trapzoidal area: 0.5*(f(x0)+f(x1))*(x1-x0)
  days_low<-seq(0,screen_interval*tiny_delta-1,by=1)/tiny_delta
  days_hi<-days_low+1/tiny_delta
  F_by_day<-0.5*(dweibull(days_low,shape=weibull_shape,scale=dwell_scale)+dweibull(days_hi,shape=weibull_shape,scale=dwell_scale))
  slip <- sum(F_by_day)*(1/(tiny_delta*screen_interval)) #day width in years
  slip
}

#use integrated weibull cumulative distribution functions
#"exact" solution to slip rate

#' @title Compute slip rate from dwell time
#' @description Compute slip rate from dwell time
#' @param dwell_model_all_df data frame with dwell time, scenario, and weibull shape
#' @param screen_interval screening interval in years
#' @param weibull_shape weibull shape parameter
#' @details 
#' dwell_model_all_df is a data frame with columns:
#' | Cancer | dwell_group | scenario | Stage | dwell | number_stage |
#' |--------|-------------|----------|-------|-------|--------------|
#'   | Anus   | A           | 1        | I     | 4     | 1            |
#'   | Anus   | A           | 1        | II    | 2     | 2            |
#'   | Anus   | A           | 1        | III   | 1     | 3            |
#'   | Anus   | A           | 1        | IV    | 1     | 4            |
#'   | Anus   | A           | 2        | I     | 3     | 1            |
#'   | Anus   | A           | 2        | II    | 1.5   | 2            |
#'   | Anus   | A           | 2        | III   | 0.75  | 3            |
#'   | Anus   | A           | 2        | IV    | 0.75  | 4            |
#'   | Anus   | A           | 3        | I     | 2     | 1            |
#'   | Anus   | A           | 3        | II    | 1     | 2            |
  
exact_slip_rate_from_dwell<-function(dwell_model_all_df,screen_interval=1,weibull_shape=1){
  # browser()
   #slip is "before clinical"
  #slip_clinical = "at stage of clinical detection"
  #assume expected is half-duration of stage of clinical detection
  #completeness in modeling
  dwell_slip_df<-dwell_model_all_df %>%
    mutate(slip = sapply(dwell,function(z){integrate_slip_rate(screen_interval,weibull_shape,z)}),
           slip_clinical=sapply(dwell*0.5,function(z){integrate_slip_rate(screen_interval,weibull_shape,z)}),
           screen_interval=screen_interval)
  
  dwell_ideal_df <- dwell_slip_df %>%
    filter(scenario==1) %>%
    mutate(dwell=10000,slip=0,slip_clinical=0,scenario=as.integer(0),screen_interval=NA)
  
  #add "scenario 0" perfect interception
  dwell_slip_df<-bind_rows(dwell_slip_df,dwell_ideal_df)
  dwell_slip_df
}
