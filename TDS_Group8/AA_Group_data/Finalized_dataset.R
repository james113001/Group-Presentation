library(forcats)
library(dplyr)
library(tidyr)
library(data.table)
library(missForest)
library(ggplot2)
library(fastDummies)
rm(list=ls())

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outputs")
data<-readRDS("./official_dataset.rds")

data<-data
#understand the values:
annot<-(readRDS("./annot.rds"))



#####################
## clean the variables which may have outliers:

continuous_data<-c("leukocyte_count",                   
                   "lymp_perc","lymp_count","mon_count","Cystatin_C",
                   "NO_poll","Z_adj_TS") #Add telomeres 

vector<-matrix(NA,nrow=length(continuous_data),ncol=2)

range_valid<-lapply(data %>% select(continuous_data[1:length(continuous_data)]),FUN= function(x){
  a <- IQR(x,na.rm=TRUE)*1.5
  vector <- mean(x,na.rm=TRUE)+c(-a,a)  
})

cts_dataset_cleanup <- data %>% select((continuous_data)) 

result <- cts_dataset_cleanup %>% mutate(
  leukocyte_count= ifelse(leukocyte_count<range_valid$leukocyte_count[1] | leukocyte_count>range_valid$leukocyte_count[2], NA, leukocyte_count ),
  lymp_perc = ifelse(lymp_perc<range_valid$lymp_perc[1] | lymp_perc>range_valid$lymp_perc[2], NA, lymp_perc ),
  NO_poll = ifelse(NO_poll<range_valid$NO_poll[1] | NO_poll>range_valid$NO_poll[2], NA, NO_poll ),
  lymp_count = ifelse(lymp_count<range_valid$lymp_count[1] | lymp_count>range_valid$lymp_count[2], NA, lymp_count ),
  Cystatin_C = ifelse(Cystatin_C<range_valid$Cystatin_C[1] | Cystatin_C>range_valid$Cystatin_C[2], NA, Cystatin_C ),
  mon_count = ifelse(mon_count<range_valid$mon_count[1] | mon_count>range_valid$mon_count[2], NA, mon_count ))

data <- cbind(result,data %>% select(-(continuous_data)))

#####################
# Filter the dataset:


df<- data %>% filter(is.na(Z_adj_TS)==FALSE) %>% 
  select(-c("Age_recruited","smoking_status","Unadj_TS","Adj_TS",
            "TS_reg_dil"))%>% 
  mutate(
  IMD = (ifelse(is.na(IMD_england),
               ifelse(is.na(IMD_scotland),IMD_wales,IMD_scotland),
               IMD_england)),
  
  Location=as.factor(ifelse(is.na(IMD_england),
                  ifelse(is.na(IMD_scotland),"Wales","Scotland"),
                  "England")),
  
  birth_weight = ifelse(birth_weight>5,birth_weight/2.205,birth_weight),
  
  distance_home_work = ifelse(distance_home_work>=1000,100,distance_home_work),
  
  dur_mod_exc = ifelse(dur_mod_exc/1440>=0.75,1440*0.75,dur_mod_exc),
  
  dur_vig_exc = ifelse(dus_vig_exc/1440>=0.75,1440*0.75,dus_vig_exc),
  
  SleepDur = ifelse(SleepDur>15,15,SleepDur),
  
  Age = 2021 - YOB,
  
  Vascular_disease=as.factor(Vascular_disease),
  
  heart_attack_bin =  as.factor(ifelse(!is.na (age_heart_attack), 1, 0)),
  
  medic_comorbidities_bin= as.factor(ifelse(!is.na (medic_comorbidities), 1, 0)),
  
  Age_diabetes_bin =  as.factor(ifelse(!is.na (Age_diabetes), 1, 0)),
  
  MenopauseAge_bin = as.factor(ifelse(!is.na (MenopauseAge), 1, 0)),
  
  Age_HBP_bin=as.factor(ifelse(!is.na (Age_HBP), 1, 0)),
  
  breadcons =ifelse(breadcons>100,100,breadcons),
  
  tea_intake =ifelse(tea_intake>30,30,tea_intake),
  
  Water_intake=ifelse(Water_intake>30,30,tea_intake),
  
  Sunburn_childhood=ifelse(Sunburn_childhood>100,100,Sunburn_childhood),
  
  Age_HBP=ifelse(Age_HBP>100,100,Age_HBP),
  
  Dia_BP=Dia_BP_auto,
  Syst_BP=Syt_BP_auto,
  
  Age_HBP=ifelse(is.na(medic_comorbidities),0,medic_comorbidities),
  
  numb_household = ifelse(numb_household>30,30,numb_household),
  
  avg_beer = ifelse(avg_beer>60,60,avg_beer),
  
  Sunburn_childhood = ifelse(is.na(Sunburn_childhood),0,Sunburn_childhood),
  
  myocardial_binary = ifelse(is.na(myocardial_date),0,1),
  
  dementia_binary = ifelse(is.na(stdate_dementia),0,1),
  
  cancer_bin=as.factor(ifelse(!is.na (cancer_date), 1, 0)),
  
  Location = as.factor(Location))

summary(df[,!sapply(df,is.factor)])

################################################################################################################################################################################################################
# Now work on factors:


#Add ETHNICITY AND VASCULAR DISEASe, cyontine thing

df_1<-df %>% mutate(
  
  Vascular_disease = recode(Vascular_disease,
                    "None of the above" = "Not disease",
                    "Prefer not to answer" = "DN/PntA",
                    "Heart attack" = "Heart attack/Angina/Stroke",
                    "Stroke" = "Heart attack/Angina/Stroke",
                    "Angina" = "Heart attack/Angina/Stroke",
                    "High blood pressure"="HBP"),
  Vascular_disease = fct_recode(Vascular_disease, NULL = "DN/PntA"),
  
  WeedFreq = recode(WeedFreq,"Prefer not to answer" = "DN/PntA",
                    "Do not know" = "DN/PntA",
                    "Less than once a month" = "Less than once a month",
                    "Once a month or more, but not every week"= "Once a month or more, but not every week",
                    "Once a week or more, but not every day"="Weekly/daily basis",
                    "Every day"="Weekly/daily basis"),
  
                  LastWeedAge = ifelse(is.na(LastWeedAge),0,LastWeedAge),
  
  souce_myocard = fct_explicit_na(souce_myocard, "Not diagnosed"),
  
  souce_myocard= recode(souce_myocard,"Self-reported only"="Self-reported only",
                        "Hospital admission"="Hospital admission/Death",
                        "Death only"="Hospital admission/Death",
                        "Not diagnosed"="Not diagnosed"),
  source_dementia = fct_explicit_na(source_dementia, "Not diagnosed"),
  
  source_dementia = recode(source_dementia,
                           "Self-reported only"="Self-reported only",
                           "Hospital admission"="Hospital admission/Death",
                           "Death only"="Hospital admission/Death",
                           "Not diagnosed"="Not diagnosed"),
  source_dementia=fct_recode(source_dementia, NULL = "Self-reported only"),
  
  alcohol_drink_status = fct_recode(alcohol_drink_status, NULL = "Prefer not to answer"),
  
  illnesses_mother = fct_explicit_na(illnesses_mother, "DN/PntA"),
  illnesses_mother= recode(illnesses_mother,
                           "None of the above (group 2)"="None of the above",
                           "None of the above (group 1)"="None of the above",
                           "Do not know (group 1)"="DN/PntA",
                           "Do not know (group 2)"="DN/PntA",
                           "Prefer not to answer (group 2)"="DN/PntA",
                           "Prefer not to answer (group 1)"="DN/PntA",
                           "DN/PntA"="DN/PntA",
                           "Heart disease"="Heart disease/Stroke",
                           "Stroke"="Heart disease/Stroke",
                           "Lung cancer"="cancer",
                           "Bowel cancer"="cancer",
                           "Breast cancer"="cancer",
                           "Chronic bronchitis/emphysema"="Other",
                           'High blood pressure'="High blood pressure",
                           "Diabetes"="Other",
                           "Alzheimer's disease/dementia"="Other",
                           "Parkinson's disease"="Other",
                           "Severe depression"="Other",
                           "Prostate cancer"="cancer",
                           "DN/PntA"="DN/PntA"),
  
  illnesses_father = fct_explicit_na(illnesses_father, "DN/PntA"),
  illnesses_father= recode(illnesses_father,
                           "None of the above (group 2)"="None of the above",
                           "None of the above (group 1)"="None of the above",
                           "Do not know (group 1)"="DN/PntA",
                           "Do not know (group 2)"="DN/PntA",
                           "Prefer not to answer (group 2)"="DN/PntA",
                           "Prefer not to answer (group 1)"="DN/PntA",
                           "DN/PntA"="DN/PntA",
                           "Heart disease"="Heart disease/Stroke",
                           "Stroke"="Heart disease/Stroke",
                           "Lung cancer"="cancer",
                           "Bowel cancer"="cancer",
                           "Breast cancer"="cancer",
                           "Chronic bronchitis/emphysema"="Other",
                           'High blood pressure'="High blood pressure",
                           "Diabetes"="Other",
                           "Alzheimer's disease/dementia"="Other",
                           "Parkinson's disease"="Other",
                           "Severe depression"="Other",
                           "Prostate cancer"="cancer",
                           "DN/PntA"="DN/PntA"),
  
  illnesses_siblings = fct_explicit_na(illnesses_siblings, "DN/PntA"),
  
  illnesses_siblings= recode(illnesses_siblings,
                           "None of the above (group 2)"="None of the above",
                           "None of the above (group 1)"="None of the above",
                           "Do not know (group 1)"="DN/PntA",
                           "Do not know (group 2)"="DN/PntA",
                           "Prefer not to answer (group 2)"="DN/PntA",
                           "Prefer not to answer (group 1)"="DN/PntA",
                           "DN/PntA"="DN/PntA",
                           "Heart disease"="Heart disease/Stroke",
                           "Stroke"="Heart disease/Stroke",
                           "Lung cancer"="cancer",
                           "Bowel cancer"="cancer",
                           "Breast cancer"="cancer",
                           "Chronic bronchitis/emphysema"="Other",
                           'High blood pressure'="High blood pressure",
                           "Diabetes"="Other",
                           "Alzheimer's disease/dementia"="Other",
                           "Parkinson's disease"="Other",
                           "Severe depression"="Other",
                           "Prostate cancer"="cancer",
                           "DN/PntA"="DN/PntA"),
  
  n_smoke_hshld= recode(n_smoke_hshld,
                      "Prefer not to answer" = "DN/PntA",
                      "No"= " No",
                      "Yes, one household member smokes"="Yes" ,
                      "Yes, more than one household member smokes"="Yes"),
  n_smoke_hshld = fct_explicit_na(n_smoke_hshld, "DN/PntA"),
  PastSmok=fct_recode(PastSmok, NULL = "DN/PntA"),

  coffee_type_intake = recode(coffee_type_intake,
                              "Prefer not to answer" = "DN/PntA",
                              "Do not know" = "DN/PntA",
                              "Decaffeinated coffee (any type)" = "Decaffeinated coffee",
                              "Instant coffee" = "Instant",
                              "Ground coffee (include espresso, filter etc)"="Ground coffee",
                              "Other type of coffee"="Other"),
  coffee_type_intake = fct_explicit_na(coffee_type_intake, "DN/PntA"),
  
  mat_smoking_birth= recode(mat_smoking_birth,
                        "Prefer not to answer" = "DN/PntA",
                        "Do not know" = "DN/PntA",
                        "No"="No",
                        "Yes"="Yes"),
#create new variable
  MenopauseAge_1 =  cut_number(MenopauseAge,n=3),
#Combine variables
  MenopauseAge =case_when(!is.na(MenopauseAge) ~ as.character(MenopauseAge_1),
                          is.na(MenopauseAge) & Sex == "Male" ~ "Male",
                          is.na(MenopauseAge) & Sex != "Male" ~ "Not yet"),
#finish
  MenopauseAge=as.factor(MenopauseAge))

print("the results are so far df_1 good")
##########################################################################################

df_2<-  df_1 %>% mutate( 
  #create new variable
  num_stillbirths_1 =  as.factor(ifelse(num_stillbirths>0,"Yes","None")),
  #Combine variables
  num_stillbirths =case_when(!is.na(num_stillbirths) ~ as.character(num_stillbirths_1),
                          is.na(num_stillbirths) & Sex == "Male" ~ "Male",
                          is.na(num_stillbirths) & Sex != "Male" ~ "None"),
  #finish
  num_stillbirths=as.factor(num_stillbirths),
  
  #create new variable
  num_miscarriages_1 =  as.factor(ifelse(num_miscarriages>0,"Yes","None")),
  #Combine variables
  num_miscarriages =case_when(!is.na(num_miscarriages) ~ as.character(num_miscarriages_1),
                             is.na(num_miscarriages) & Sex == "Male" ~ "Male",
                             is.na(num_miscarriages) & Sex != "Male" ~ "None"),
  #finish
  num_miscarriages=as.factor(num_miscarriages),
  
  contra_pill = case_when(Sex == "Male"~ "Male",
                          Sex != "Male"~ as.character(contra_pill)),
  
  contra_pill = recode(contra_pill,
    "None of the above"= "DN/PntA",
    "Prefer not to answer"="DN/PntA",
    "No"="No",
    "Yes" = "Yes"),
    
  hormone_rep_thera = case_when(Sex == "Male"~ "Male",
                                Sex != "Male"~ as.character(hormone_rep_thera)),
  
  hormone_rep_thera = recode(hormone_rep_thera,
                             "None of the above"= "DN/PntA" ,
                             "Prefer not to answer"="DN/PntA", 
                             "No"="No",
                             "Yes"="Yes"),
  
  menopause_bi =case_when(Sex == "Male"~ "Male",
                          Sex != "Male"~ as.character(menopause_bi)),
  
  menopause_bi = recode(menopause_bi,
    "None of the above"= "DN/PntA" ,
    "Prefer not to answer"="DN/PntA",
    "Not sure - had a hysterectomy" = "Not sure",
    "Not sure - other reason" = "Not sure",
    "Male" = "Male"),
  
  menopause_bi=fct_recode(menopause_bi, NULL = "DN/PntA"),
  
  presc_medicine= recode(presc_medicine,
                         "Do not know"= "DN/PntA" ,
                         "Prefer not to answer"="DN/PntA",
                         "No"="No",
                         "Yes - you will be asked about this later by an interviewer"="Yes"))

##########################################################################################

df_3<- df_2 %>%mutate(  accomodation_lived = recode(accomodation_lived,
     "None of the above"= "Other" ,
     "Prefer not to answer"="Other" ,
     "A house or bungalow"="House/Bungalow" ,
     "A flat, maisonette or apartment"="Flat, maisonette or apartment" ,
     "Mobile or temporary structure (i.e. caravan)"="Other" ,
     "Sheltered accommodation"="Other" ,
     "Care home"="Other" 
  ),
  
  own_accomodation = recode(own_accomodation,
    "Own outright (by you or someone in your household)"="Own" ,
    "Own with a mortgage"="Own",
    "Pay part rent and part mortgage (shared ownership)" = "Own",
    "Rent - from local authority, local council, housing association"="Rent",
    "Rent - from private landlord or letting agency"="Rent" ,
    "Live in accommodation rent free"="Other",
    "Prefer not to answer"="DN/PntA",
    "(Other)"="Other",
    "None of the above"="Other"),
  own_accomodation=fct_recode(own_accomodation, NULL = "DN/PntA"),
  
  avg_household_income = recode(avg_household_income,
    "Prefer not to answer"="DN/PntA",
    "Do not know"="DN/PntA",
    "Less than 18,000" = "Less than 18,000",
    "18,000 to 30,999"= "18,000 to 30,999",
    "31,000 to 51,999" = "31,000 to 51,999",
    "52,000 to 100,000" = "52,000 to 100,000",
    "Greater than 100,000" ="Over 100,000"),
  
  walk_pace = recode(walk_pace,
   "Prefer not to answer" = "DN/PntA",
    "Do not know"= "DN/PntA",
    "None of the above"="DN/PntA",
    "Slow pace" ="Slow pace",
    "Steady average pace"= "Steady average pace",
    "Brisk pace" = "Brisk pace"),
  walk_pace=fct_recode(walk_pace, NULL = "DN/PntA"),
  
  Insomnia = recode(Insomnia,
                    "Prefer not to answer"="DN/PntA",
                    "Never/rarely"="Never/rarely",
                    "Sometimes"="Sometimes",
                    "Usually"="Usually"),
  Insomnia=fct_recode(Insomnia, NULL = "DN/PntA"),
  
  narcolepsy = recode(narcolepsy,
                     "Prefer not to answer"="DN/PntA",
                     "Do not know"="DN/PntA",
                     "Never/rarely"= "Never/rarely",
                     "Sometimes"= "Sometimes",
                     "Often"= "Often/Always",
                     "All of the time"= "Often/Always"),
  narcolepsy=fct_recode(narcolepsy, NULL = "DN/PntA"),
  
  CurrSmok= recode(CurrSmok,
                   "Yes, on most or all days"="Yes",
                  "Only occasionally" = "Yes" ,
                   "No" = "No",
                   "Prefer not to answer"= "DN/PntA"),
  CurrSmok=fct_recode(CurrSmok, NULL = "DN/PntA"),
  
  PastSmok= recode(PastSmok,
                   "Prefer not to answer"= "DN/PntA",
                   "Smoked on most or all days"="Smoked on most or all days",
                   "Smoked occasionally"="Smoked occasionally",
                   "Just tried once or twice"="Never",
                   "I have never smoked"="Never"),
  PastSmok=fct_recode(PastSmok, NULL = "DN/PntA"),
  
  oily_fish_intake= recode(oily_fish_intake,
                           "Prefer not to answer"= "DN/PntA",
                           "Do not know"="DN/PntA",
                           "Never"="Never",
                           "Less than once a week"="Less or once a week",
                           "Once a week"= "Less or once a week",
                           "2-4 times a week"="2-4 times a week",
                           "5-6 times a week"="5-6 times a week/Daily",
                           "Once or more daily"="5-6 times a week/Daily"),
  oily_fish_intake=fct_recode(oily_fish_intake, NULL = "DN/PntA"),
  
  non_oily_fish_intake= recode(non_oily_fish_intake,
                               "Prefer not to answer"= "DN/PntA",
                               "Do not know"="DN/PntA",
                               "Never"="Never",
                               "Less than once a week"="Less or once a week",
                               "Once a week"= "Less or once a week",
                               "2-4 times a week"="2-4 times a week",
                               "5-6 times a week"="5-6 times a week/Daily",
                               "Once or more daily"="5-6 times a week/Daily"),
  non_oily_fish_intake=fct_recode(non_oily_fish_intake, NULL = "DN/PntA"),
  
  ProcMeat= recode(ProcMeat,
                   "Prefer not to answer"= "DN/PntA",
                   "Do not know"="DN/PntA",
                   "Never"="Never",
                   "Less than once a week"="Less or once a week",
                   "Once a week"= "Less or once a week",
                   "2-4 times a week"="2-4 times a week",
                   "5-6 times a week"="5-6 times a week/Daily",
                   "Once or more daily"="5-6 times a week/Daily"),
  ProcMeat=fct_recode(ProcMeat, NULL = "DN/PntA"),
  
  cheese_cons= recode(cheese_cons,
                      "Prefer not to answer"= "DN/PntA",
                      "Do not know"="DN/PntA",
                      "Never"="Never",
                      "Less than once a week"="Less or once a week",
                      "Once a week"= "Less or once a week",
                      "2-4 times a week"="2-4 times a week",
                      "5-6 times a week"="5-6 times a week/Daily",
                      "Once or more daily"="5-6 times a week/Daily"),
  cheese_cons=fct_recode(cheese_cons, NULL = "DN/PntA"),
  
  salt_in_food = recode(salt_in_food,
                        "Prefer not to answer"= "DN/PntA",
                        "Nerver/rarely"="Nerver/rarely",
                        "Sometimes"="Sometimes",
                        "Usually"="Usually",
                        "Always"="Always"),
  salt_in_food=fct_recode(salt_in_food, NULL = "DN/PntA"),
  
  Alcohol_intake_freq=fct_recode(Alcohol_intake_freq, NULL = "Prefer not to answer"),
  
  depressed_last_2wks= recode(dep_last_2_weeks,
                           "Prefer not to answer"="DN/PntA",
                           "Do not know"="DN/PntA",
                           "Not at all"="Not at all",
                           "Several days"="Several days",
                           "More than half the days"="More than half/nearly always",
                           "Nearly every day"="More than half/nearly always" ),
  dep_last_2_weeks=fct_recode(dep_last_2_weeks, NULL = "DN/PntA"),
  
    contra_pill=recode(contra_pill,
                     "Prefer not to answer"="DN/PntA",
                     "Do not know"="DN/PntA",
                     "No"= "No",
                     "Yes"="Yes"),
  contra_pill=fct_recode(contra_pill, NULL = "DN/PntA"),

  long_standing_illness= recode(long_standing_illness,
                  "Prefer not to answer"="DN/PntA",
                  "Do not know"="DN/PntA",
                  "No"= "No",
                  "Yes"="Yes"),
  long_standing_illness=fct_recode(long_standing_illness, NULL = "DN/PntA"),
  
  diabetes_doctor= recode(diabetes_doctor,
                          "Prefer not to answer"="DN/PntA",
                          "Do not know"="DN/PntA",
                          "No"= "No",
                          "Yes"="Yes"),
  diabetes_doctor=fct_recode(diabetes_doctor, NULL = "DN/PntA"),
  
  hormone_rep_thera= recode(hormone_rep_thera,
                            "Prefer not to answer"="DN/PntA",
                            "Do not know"="DN/PntA",
                            "No"= "No",
                            "Yes"="Yes"),
  hormone_rep_thera=fct_recode(hormone_rep_thera, NULL = "DN/PntA"),
  
  
  cancer_doctor = recode(cancer_doctor,
                         "Prefer not to answer"="DN/PntA",
                         "Do not know"="DN/PntA",
                         "No"= "No",
                         "Yes - you will be asked about this later by an interviewer"="Yes"),
  
  hormone_rep_thera= as.factor(hormone_rep_thera),
  menopause_bi= as.factor(menopause_bi),
  contra_pill= as.factor(contra_pill),
  Ethnicity= recode(Ethnicity,
                    "Prefer not to answer"="DN/PntA",
                    "Do not know"="DN/PntA",
                    "White"="White/White British",
                    "British"="White/White British",
                    "Irish"="White/White British",
                    "Any other white background"="White/White British",
                    "Mixed"="Mixed",
                    "White and Black Caribbean"="Mixed",
                    "White and Black African"="Mixed",
                    "White and Asian"="Mixed",
                    "Any other mixed background"="Mixed",
                    "Asian or Asian British"="Asian/Asian British",
                    "Indian"="Asian/Asian British",
                    "Pakistani"="Asian/Asian British",
                    "Bangladeshi"="Asian/Asian British",
                    "Any other Asian background"="Asian/Asian British",
                    "Black or Black British"="Black/Black British",
                    "Caribbean" = "Black/Black British",
                    "African"="Black/Black British",
                    "Any other Black background"="Black/Black British",
                    "Chinese"="Chinese",
                    "Other ethnic group"="Other"),
  Ethnicity=fct_recode(Ethnicity, NULL = "DN/PntA")
  ) %>% select(-c("YOB","IMD_england","IMD_england","IMD_wales",
                  "IMD_scotland","cancer_date","stdate_dementia",
                  "myocardial_date","dep_last_2_weeks","MenopauseAge_1","num_stillbirths_1","num_miscarriages_1",
                  "IMD_town","lymp_perc")) %>%
  mutate(
    accomodation_lived=relevel(accomodation_lived,ref="House/Bungalow"),
    own_accomodation=relevel(own_accomodation,ref="Own"),
    narcolepsy=relevel(narcolepsy,ref="Never/rarely"),
    diabetes_doctor=relevel(diabetes_doctor,ref="No"),
    cancer_doctor=relevel(cancer_doctor,ref="No"),
    presc_medicine=relevel(presc_medicine,ref="No"),
    menopause_bi=relevel(menopause_bi,ref="No"),
    contra_pill=relevel(contra_pill,ref="No"),
    hormone_rep_thera=relevel(hormone_rep_thera,ref="No"),
    num_stillbirths=relevel(num_stillbirths,ref="Male"),
    souce_myocard=relevel(souce_myocard,ref="Not diagnosed"),
    source_dementia=relevel(source_dementia,ref="Not diagnosed"),
    depressed_last_2wks=relevel(depressed_last_2wks,ref="Not at all"))




















na_perc<-round(100*sapply(df_3, function(x) sum(length(which(is.na(x)))))/length(df_3$Sex),3)


optimal_deletion<-matrix(NA,nrow=length(seq(5,10,by=0.5)),ncol=2)
j=1

for ( i in seq(5,10.5,by=0.5)){
  optimal_deletion[j,1]<-length(colnames(df_3))-length(names(na_perc[na_perc>i]))
  not_selected<-names(na_perc[na_perc>i])
  optimal_deletion[j,2] <- nrow(df_3 %>% select(-not_selected) %>% na.omit())/1000
  j=j+1
}

optimal_deletion  #Optimal value to remove is at 6% missingness === 7% level


na_perc_colnames_6<-names(na_perc[na_perc>6])
na_perc_colnames_10.5<-names(na_perc[na_perc>10.5])  

library(tidyr)

dum<-dummy_cols(df_final[,sapply(df_final,is.factor)] , select_columns = colnames(df_final[,sapply(df_final,is.factor)]))%>%
  select(-c(colnames(df_final[,sapply(df_final,is.factor)])))

dummy_columns<-cbind(dum,df_final[,!sapply(df_final,is.factor)])
saveRDS(dummy_columns, "../AA_Group_data/DATASET/ukbiobank_dummy_cols.rds")

df_final <- df_3 %>% select(-(na_perc_colnames_6)) %>% na.omit()
saveRDS(df_final, "../AA_Group_data/DATASET/data_ukbiobank.rds")

df_trees <- df_3 %>% select(-(na_perc_colnames_10.5))
saveRDS(df_trees, "../AA_Group_data/DATASET/trees_dataset.rds")

#######

#Variable groups:

lifestyle<-c("num_days_walk_10mins","num_days_mod_exc","num_says_vig_exc","walk_pace","SleepDur","Insomnia","narcolepsy",
             "CurrSmok","PastSmok","n_smoke_hshld","tobacco_exp_hshld","Location")
lifestyle<-data.frame(variables=lifestyle,type=rep("lifestyle",length(lifestyle)))

dietary_lifestyle<-c("cooked_veggie_intake","raw_veggie_intake","fresh_fruit_intake",
             "dried_fruit_intake","oily_fish_intake","non_oily_fish_intake","ProcMeat","cheese_cons","breadcons",
             "salt_in_food","tea_intake","coffee_type_intake","Water_intake","Alcohol_intake_freq","vit_suplements")
dietary_lifestyle<-data.frame(variables=dietary_lifestyle,type=rep("dietary_lifestyle",length(dietary_lifestyle)))

socio_economic<-c("IMD_town","accomodation_lived","own_accomodation","length_current_address","avg_household_income",
                  "disability_allowance","LastWeedAge","numb_household","IMD")
socio_economic<-data.frame(variables=socio_economic,type=rep("socio_economic",length(socio_economic)))

medical<-c("leukocyte_count","NO_poll","num_cancers","Sunburn_childhood","mat_smoking_birth","long_standing_illness",
           "diabetes_doctor","cancer_doctor","presc_medicine","menopause_bi","contra_pill","hormone_rep_thera","Age_HBP",
           "forces_vital_capacity","MenopauseAge","num_stillbirths","num_miscarriages","illnesses_father","illnesses_mother",
           "illnesses_siblings","NoEDWS","common_disease_diag","NO2_poll","PM10_poll","souce_myocard",
           "source_dementia","Syt_BP_auto","Dia_BP_auto","heart_attack_bin","medic_comorbidities_bin","Age_diabetes_bin",
           "MenopauseAge_bin","Age_HBP_bin","myocardial_binary","dementia_binary","depressed_last_2wks","Cystatin_C","Vascular_disease")

medical<-data.frame(variables=medical,type=rep("medical",length(medical)))
  
personal<-c("Sex","Age_recruited","BMI","Z_adj_TS","Height","Age","Ethnicity")
personal<-data.frame(variables=personal,type=rep("personal",length(personal)))

subgroup_organization<-left_join(data.frame(variables=colnames(df_final)),rbind(lifestyle,dietary_lifestyle,socio_economic,medical,personal))
saveRDS(subgroup_organization, "../AA_Group_data/DATASET/ukbiobank_subgrouped.rds")


dum<-dummy_cols(df_final[,sapply(df_final,is.factor)] , select_columns = colnames(df_final[,sapply(df_final,is.factor)]))%>%
  select(-c(colnames(df_final[,sapply(df_final,is.factor)])))

print("created the dummy file")

dummy_columns<-cbind(dum,df_final[,!sapply(df_final,is.factor)])
saveRDS(dummy_columns, "../AA_Group_data/DATASET/ukbiobank_dummy_cols.rds")

print("saved the dummy file")
