#join SEER18 data
library(tidyverse)
library(readxl)

#date_code<-"20200520"

SEER_spreadsheet_path<-"data/SEER_Draw_Interception_Manuscript_20200728.xlsx"

#should make sure to output the raw data tables from the spreadsheet for future use
# ir: incidence rate
# aair: age-adjusted incidence rate
stage_aair_data<-read_excel(SEER_spreadsheet_path,
                            sheet="Incidence",
                            range="A2:E126",
                            col_names=c("SEER_Draw","Stage","IR","Count","Population"))

stage_css_data<-read_excel(SEER_spreadsheet_path,
                           sheet="CSS",
                           range="A4:E2503",
                           col_names=c("SEER_Draw","Stage","TIME","N","CSS"),
                           guess_max=5000)

#filter out "ERROR" codes in SEER results and turn them into NA values
#and return CSS to double type
# SEER_Draw Stage TIME      N   CSS
# 1 Anus      I     3 mo   1668 0.990
# 2 Anus      I     6 mo   1668 0.984
# 3 Anus      I     9 mo   1668 0.98 
# 4 Anus      I     12 mo  1668 0.973
# 5 Anus      I     15 mo  1668 0.966
# 6 Anus      I     18 mo  1668 0.958
# 7 Anus      I     21 mo  1668 0.954
# 8 Anus      I     24 mo  1668 0.946
# 9 Anus      I     27 mo  1668 0.942
# 10 Anus     I     30 mo  1668 0.939
stage_css_data<-stage_css_data %>%
  mutate(CSS=case_when(grepl("ERROR",CSS)~NA_character_,
                       TRUE ~ CSS)) %>%
  type_convert()

# SEER_Draw Stage           TIME      N   CSS Survival
# 1 Anus      I               60 mo  1668 0.909    0.909
# 2 Anus      II              60 mo  2588 0.821    0.821
# 3 Anus      III             60 mo  2384 0.672    0.672
# 4 Anus      IV              60 mo   747 0.235    0.235
# 5 Anus      Unknown/missing 60 mo  1997 0.778    0.778
# 6 Bladder   I               60 mo 19008 0.862    0.862
# 7 Bladder   II              60 mo  8871 0.618    0.618
# 8 Bladder   III             60 mo  3554 0.495    0.495
# 9 Bladder   IV              60 mo  6716 0.163    0.163
# 10 Bladder   Unknown/missing 60 mo  2901 0.660    0.660
stage_css_filtered_data<-stage_css_data %>%
  mutate(Survival=CSS) %>%
  filter(TIME=="60 mo")

# SEER_Draw Stage              IR Survival
# 1 Anus      I                0.93    0.909
# 2 Anus      II               1.41    0.821
# 3 Anus      III              1.25    0.672
# 4 Anus      IV               0.39    0.235
# 5 Anus      Unknown/missing  1.15    0.778
# 6 Bladder   I               11.1     0.862
# 7 Bladder   II               5.05    0.618
# 8 Bladder   III              2.02    0.495
# 9 Bladder   IV               3.73    0.163
# 10 Bladder   Unknown/missing  1.95    0.660
stage_joint_filtered<-stage_aair_data %>% 
  left_join(stage_css_filtered_data) %>% 
  select(SEER_Draw,Stage,IR,Survival)

#have to deal with these cancers specially
hard_cancers<-c("Lymphoid Leukemia","Myeloid Neoplasm","Plasma Cell Neoplasm","[OTHER]")

#okay, deal with typical cases where stage exists, and unknown/missing can be imputed sensibly
# SEER_Draw Stage              IR Survival
# 1 Anus      I                0.93    0.909
# 2 Anus      II               1.41    0.821
# 3 Anus      III              1.25    0.672
# 4 Anus      IV               0.39    0.235
# 5 Anus      Unknown/missing  1.15    0.778
# 6 Bladder   I               11.1     0.862
# 7 Bladder   II               5.05    0.618
# 8 Bladder   III              2.02    0.495
# 9 Bladder   IV               3.73    0.163
# 10 Bladder   Unknown/missing  1.95    0.660
limited_joint_filtered<-stage_joint_filtered %>%
  filter(!(SEER_Draw %in% hard_cancers)) 

#impute these
# UR: unknown rate
# SEER_Draw          UR
# 1 Anus             1.15
# 2 Bladder          1.95
# 3 Breast           7.46
# 4 Cervix           0.58
# 5 Colon/Rectum     9.68
# 6 Esophagus        1.79
# 7 Gallbladder      0.88
# 8 Head and Neck    6.11
# 9 Kidney           2.26
# 10 Liver/Bile-duct  5.11
unknown_joint_filtered<-limited_joint_filtered %>%
  filter(Stage=="Unknown/missing") %>%
  group_by(SEER_Draw) %>%
  summarize(UR=sum(IR,na.rm=TRUE)) %>%
  ungroup()

tStage=c("I","II","III","IV","Unknown/missing")

# URX: unknown rate by stage
# SEER_Draw Stage     IR Survival
# 1 Anus      I      1.20     0.909
# 2 Anus      II     1.82     0.821
# 3 Anus      III    1.61     0.672
# 4 Anus      IV     0.503    0.235
# 5 Bladder   I     12.1      0.862
# 6 Bladder   II     5.50     0.618
# 7 Bladder   III    2.20     0.495
# 8 Bladder   IV     4.06     0.163
# 9 Breast    I     90.0      0.982
# 10 Breast    II    56.1      0.929
imputed_joint_filtered<-limited_joint_filtered %>%
  filter(Stage %in% tStage[1:4]) %>%
  left_join(unknown_joint_filtered) %>%
  group_by(SEER_Draw) %>%
  # attributes unknown/missing to each known stage by proportion to incidence
  mutate(URX=UR*IR/sum(IR,na.rm=TRUE)) %>%
  mutate(URX=replace_na(URX,0.0)) %>%
  ungroup() %>%
  mutate(IR=IR+URX) %>%
  select(-UR,-URX)

#unstaged and expected not to be staged
#need to up-impute "staged" to "notstaged" for lymphoid leukemia
#because we don't have by-stage sensitivities that are relevant to those entries in SEER
#rate is relatively small

# treat Lympoid Leukemia, Myeloid Neoplasm, Plasma Cell Neoplasm as unstaged and IR is the total incidence from all stages (if available)
# SEER_Draw            Stage        IR Survival
# 1 Lymphoid Leukemia    NotStaged  17.4    0.853
# 2 Myeloid Neoplasm     NotStaged  14.6    0.431
# 3 Plasma Cell Neoplasm NotStaged  17.9    0.593
unstaged_joint_filtered<-stage_joint_filtered %>%
  filter(SEER_Draw %in% hard_cancers[1:3]) %>%
  group_by(SEER_Draw) %>%
  summarize(
    IR=sum(IR,na.rm=TRUE),
    Survival=Survival[Stage=="Unknown/missing"]) %>%
  ungroup() %>%
  mutate(Stage="NotStaged") %>%
  select(SEER_Draw,Stage,IR,Survival)

#other - heterogenous group so cannot impute unstaged to staged
#also no sensible group-level sensitivity
#but do need incidence and survival data
#   SEER_Draw Stage        IR Survival
# 1 [OTHER]   I          4.03    0.857
# 2 [OTHER]   II         2.14    0.705
# 3 [OTHER]   III        2.46    0.540
# 4 [OTHER]   IV         2.33    0.199
# 5 [OTHER]   NotStaged 59.7     0.429
other_joint_filtered<-stage_joint_filtered %>%
  filter(SEER_Draw=="[OTHER]") %>%
  mutate(Stage=case_when(Stage!="Unknown/missing" ~ Stage,
                         TRUE ~ "NotStaged"))

# SEER_Draw Stage     IR Survival
# 1 Anus      I      1.20     0.909
# 2 Anus      II     1.82     0.821
# 3 Anus      III    1.61     0.672
# 4 Anus      IV     0.503    0.235
# 5 Bladder   I     12.1      0.862
# 6 Bladder   II     5.50     0.618
# 7 Bladder   III    2.20     0.495
# 8 Bladder   IV     4.06     0.163
# 9 Breast    I     90.0      0.982
# 10 Breast   II    56.1     0.929
total_joint_filtered<-bind_rows(imputed_joint_filtered,unstaged_joint_filtered,other_joint_filtered)

#output imputed SEER data
write_tsv(total_joint_filtered,sprintf("generated_data/%s_total_seer_draw.tsv",date_code))


