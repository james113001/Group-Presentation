rm(list=ls())
library(tidyr)
library(dplyr)
library(data.table)

crit <- data.frame(fread("/rds/general/project/hda_21-22/live/TDS/General/Data/hesin_critical.txt")) %>% mutate(eid=as.character(eid))

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outcome_definitions_COVID")
covid_deaths<-readRDS("./Outputs/output_final.rds")%>% mutate(eid=as.character(eid))
covid_tests<-readRDS("./Outputs/covid_tested.rds")%>% mutate(eid=as.character(eid))



################################################################################################
## Initial information data

deathEID <- covid_deaths[,1]
covidEID <- covid_tests[,1]

length(unique(crit[["eid"]])) #number people in ICU
length(unique(covid_tests[["eid"]])) #168360 num people in tests

#mutate eid for later:

covid_deaths<-covid_deaths %>% mutate(eid=as.character(eid))
covid_deaths$eid<-as.integer(rownames(covid_deaths))  

covid_crit<-merge(crit, covid_tests, by = "eid")%>% mutate(
  ccstartdate=as.Date(ccstartdate, format="%d/%m/%Y"),
  specdate=as.Date(specdate, format="%d/%m/%Y"),
  ccdisdate=as.Date(ccdisdate, format="%d/%m/%Y"))

crit<-crit %>% mutate(
  ccstartdate=as.Date(ccstartdate, format="%d/%m/%Y"),
  specdate=as.Date(specdate, format="%d/%m/%Y"),
  ccdisdate=as.Date(ccdisdate, format="%d/%m/%Y"))
#n=64167 encounters

length(unique(covid_crit[["eid"]])) #8958 pt

##ccstartdate filter to Jan 31, 2020 (first COVID case in UK)
#covid_crit<-covid_crit[covid_crit$ccstartdate >= "2020-01-31"] #8274 pt, n=58699 encounters

################################################################################################

## filter the filter the first covid case
covid_crit<- covid_crit %>% filter(ccdisdate >= "2020-01-31") #8308 pt, n=58885 encounters

##          ccadmittype, filter for only 1 and 2 (emergency/unplanned):
#covid_emergency<-covid_crit[covid_crit$ccadmitype == 1 | covid_crit$ccadmitype == 2] #n=23875 encounters
#length(unique(covid_emergency[["eid"]]))#3297 pt
#covid_emergency_inf<-covid_emergency %>% filter(result == 1) #n=2320 encounters 
#length(unique(covid_emergency_inf[["eid"]])) #743pt too sparse!


##ccadmisorc, ccsorcloc?, ccdisrdydate? (redundant w. ccdisdate), ccdisdest (redundant with ccdisloc) delete
##dermsup delete
##bressupdays (basic resp support) aressupdays(advanced)
##bcardsupdays/acardsupdays (cardiac support) 
##rensupdays (renal), livsup(liver)
##neurosup delete? ppl w neuro diseases suscept. to cov infection
##GIsup (GI), orgsupmax (max # of organ sup needed in a day)
##unitbedconfig? (redundant with organ support var.) delete



################################################################################################
# Check how many if no filter by admission type:

covid_emergency_noadmitt<- covid_crit %>% filter(result == 1)
length(unique(covid_emergency_noadmitt[["eid"]])) #doubled the sample size!


#Alternatively, for each ICU, if not repeated double up and substract dates to see date difference between test/admission.

variables_discard<- c("unitbedconfig","ccunitfun","cclev3days","cclev2days","dsource","source","arr_index","ins_index","ccapcrel")

covid_emergency_noadmitt<- covid_crit %>% filter(result == 1) %>% 
  select(-(variables_discard)) %>%  # remove the variables stated above, not useful
  mutate(test_to_icu = abs(as.integer(specdate)-as.integer(ccstartdate)),
         release_to_positive =as.integer(specdate)-as.integer(ccdisdate),
         Loc = as.factor(Loc),
         icu_stay=as.factor(ifelse(is.na(ccstartdate),0,1)),
         bsc_res_sup =as.factor(ifelse(bressupdays>0,1,0)) , # 1 = Yes
         adv_res_sup= as.factor(ifelse(aressupdays>0,1,0)) , # 1 = Yes
         other_organ_support=as.factor(ifelse(bcardsupdays>0 | # if any of these variables is positive mark it
                                              acardsupdays>0| 
                                              rensupdays>0 |
                                              neurosupdays>0| 
                                              gisupdays>0| 
                                              dermsupdays>0| 
                                              liversupdays>0,1,0))) %>% 
  filter(test_to_icu<=15 & release_to_positive < 2 ) %>% select(-c(site,ccdisrdydate))

# Total number of  people who have gotten COVID after/before 15 days being in ICU, 
# but not so close that they might have caught it there and being released.

length(unique(covid_emergency_noadmitt[["eid"]])) #end up with 424 individuals in ICU only.

##################################################################
# Final cut for variables:


remove_vars<-c("ccstartdate","ccadmitype","ccadmisorc","ccdisstat","ccdisdest","ccdisloc",
               "Loc","origin","pattype","specdate","release_to_positive","test_to_icu","ccsorcloc","ccdisdate")

#THIS LINE OF CODE CUTS THE SMALLEST TEST-TO-ICU DATE, so cut repeats of eid's:

single_eid_icu_covd<-covid_emergency_noadmitt %>% group_by(eid) %>% 
  slice(which.min(test_to_icu)) %>% ungroup(eid) %>% select(-remove_vars)
  
  # combine ICU data and deaths, but use deaths data as main cases

removal_deaths<-c("date_recr","date_diagnosis","date_death","prevalent_case","incident_case","time_to_diagnosis")

deaths_cleaned <-covid_deaths %>% 
  filter(!is.na(date_death),!is.na(date_diagnosis)) %>% 
  mutate(
    death=as.factor(ifelse(is.na(date_death),0,1)),
    eid = as.character(eid)) %>% 
  select(-removal_deaths)

deaths_icu_covid<-left_join(single_eid_icu_covd,deaths_cleaned %>% 
                            select(eid,death),by="eid") %>% 
                            mutate(death = ifelse(is.na(death),0,1)) %>%
                            filter(result==1)

#Left join deaths data on ICU, to see if patients died

final_database_icu_death<-data.frame(deaths_icu_covid)
rownames(final_database_icu_death)<-deaths_icu_covid$eid

saveRDS(final_database_icu_death,"../AA_Group_data/DATASET/deaths_icu_covid.rds")

eid_positives<-unique(covid_tests %>% filter(result==1) %>% select(eid))
saveRDS(eid_positives,"../AA_Group_data/DATASET/test_covid.rds")

##################################################################
# Same procedure but with NON COVID PATIENTS


non_covid_emergency <- crit %>%
  select(-(variables_discard)) %>%  # remove the variables stated above, not useful
  mutate(Loc = as.factor(Loc),
         icu_stay=as.factor(ifelse(is.na(ccstartdate),0,1)),
         bsc_res_sup =as.factor(ifelse(bressupdays>0,1,0)) , # 1 = Yes
         adv_res_sup= as.factor(ifelse(aressupdays>0,1,0)) , # 1 = Yes
         other_organ_support=as.factor(ifelse(bcardsupdays>0 | # if any of these variables is positive mark it
                                                acardsupdays>0| 
                                                rensupdays>0 |
                                                neurosupdays>0| 
                                                gisupdays>0| 
                                                dermsupdays>0| 
                                                liversupdays>0,1,0))) %>% 
  filter(!(eid%in%covid_emergency_noadmitt$eid)) %>% select(-c(site,ccdisrdydate)) %>% select(-remove_vars)

saveRDS(non_covid_emergency,"../AA_Group_data/DATASET/non_covid_icu.rds")

deaths_icu_non_covid<-deaths_cleaned %>%
        filter(case!=1)%>% 
        select(eid,death) %>% 
        mutate(death = ifelse(is.na(death),0,1)) 

saveRDS(deaths_icu_non_covid,"../AA_Group_data/DATASET/non_covid_deaths.rds")
