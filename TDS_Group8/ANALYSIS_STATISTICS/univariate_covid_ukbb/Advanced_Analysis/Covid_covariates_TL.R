rm(list=ls())
library(tidyverse)
library(glmnet)
library(focus)
library(withr)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(data.table)
################ REAL DEATHS!!
setwd("/rds/general/project/hda_21-22/live/TDS/General/Data")
deathcause<-fread("death_cause.txt")
deathdates<-fread("death.txt")

deathdeets<-deathcause %>%
  filter(cause_icd10=="U071")%>%
  mutate(eid=as.character(eid)) %>% distinct() %>% select(eid)

################


setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis")

#Get the dataset with the extra varables, but more constrained:

ukbiobank_dataset<-readRDS("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data/DATASET/trees_dataset.rds")
colnames(ukbiobank_dataset)
#Add:: Ethnicity and cycstatin, what happened with Alcohol_drinker_status

Covid_covariates<-ukbiobank_dataset%>% 
  select("eid","Age","Ethnicity","Sex","Diabetes","High_blood_pressure","Vascular_disease","Average_household_income","leukocyte_count","Myocardial_disease","Smoking_status",
         "Alcohol_drinker_status","Leukocyte_telomere_length")%>%
  na.omit()

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")


dff1<-readRDS("DATASET/deaths_icu_covid.rds")%>% 
  mutate(eid = as.character(eid),
         case=as.factor(case),
         icu_stay=as.factor(icu_stay),
         bsc_res_sup=recode(as.factor(bsc_res_sup),
                            "0"="0",
                            "1"="1",
                            "2"="1"),
         adv_res_sup=recode(as.factor(adv_res_sup),
                            "0"="0",
                            "1"="1",
                            "2"="1"),
         other_organ_support=recode(as.factor(other_organ_support),
                                    "0"="0",
                                    "1"="1",
                                    "2"="1"))


dff2<-readRDS("DATASET/test_covid.rds") %>% 
  mutate(results=1, eid = as.character(eid))

non_covid_icu<-readRDS("DATASET/non_covid_icu.rds")%>% select(-liversupdays) %>%
  mutate(eid = as.character(eid),
         icu_stay=as.factor(icu_stay),
         bsc_res_sup=recode(as.factor(bsc_res_sup),
                            "0"="0",
                            "1"="1",
                            "2"="1"),
         adv_res_sup=recode(as.factor(adv_res_sup),
                            "0"="0",
                            "1"="1",
                            "2"="1"),
         other_organ_support=recode(as.factor(other_organ_support),
                                    "0"="0",
                                    "1"="1",
                                    "2"="1")) %>% mutate(
                                      bressupdays=ifelse(is.na(bressupdays),0,bressupdays),
                                      aressupdays=ifelse(is.na(aressupdays),0,aressupdays),
                                      bcardsupdays=ifelse(is.na(bcardsupdays),0,bcardsupdays),
                                      acardsupdays=ifelse(is.na(acardsupdays),0,acardsupdays),
                                      rensupdays=ifelse(is.na(rensupdays),0,rensupdays),
                                      neurosupdays=ifelse(is.na(neurosupdays),0,neurosupdays),
                                      gisupdays=ifelse(is.na(gisupdays),0,gisupdays),
                                      dermsupdays=ifelse(is.na(dermsupdays),0,dermsupdays),
                                      orgsupmax=ifelse(is.na(orgsupmax),0,orgsupmax),
                                      icu_stay=fct_explicit_na(icu_stay, "0"),
                                      bsc_res_sup=fct_explicit_na(bsc_res_sup, "0"),
                                      adv_res_sup=fct_explicit_na(adv_res_sup, "0"),
                                      other_organ_support=fct_explicit_na(other_organ_support, "0"))

non_covid_deaths<-readRDS("DATASET/non_covid_deaths.rds") %>% 
  mutate(eid = as.character(eid))



icu_death_ukbb<-left_join(left_join(Covid_covariates,dff1,by="eid") %>% mutate(
  death = as.factor(ifelse(is.na(death),0,death)),
  bressupdays=ifelse(is.na(bressupdays),0,bressupdays),
  aressupdays=ifelse(is.na(aressupdays),0,aressupdays),
  bcardsupdays=ifelse(is.na(bcardsupdays),0,bcardsupdays),
  acardsupdays=ifelse(is.na(acardsupdays),0,acardsupdays),
  rensupdays=ifelse(is.na(rensupdays),0,rensupdays),
  neurosupdays=ifelse(is.na(neurosupdays),0,neurosupdays),
  gisupdays=ifelse(is.na(gisupdays),0,gisupdays),
  dermsupdays=ifelse(is.na(dermsupdays),0,dermsupdays),
  orgsupmax=ifelse(is.na(orgsupmax),0,orgsupmax),
  icu_stay=fct_explicit_na(icu_stay, "0"),
  bsc_res_sup=fct_explicit_na(bsc_res_sup, "0"),
  adv_res_sup=fct_explicit_na(adv_res_sup, "0"),
  other_organ_support=fct_explicit_na(other_organ_support, "0"),
  case=fct_explicit_na(as.factor(case), "0")),
  dff2,by="eid") %>% select(-case)



# Matching on: Age, Sex, BMI, Ethnicity?
##############################################################################################################################

ICU_total<-rbind(non_covid_icu %>% mutate(case=as.factor(0),death=as.factor(0)),dff1)

ICU_ukbiobank<-left_join(left_join(Covid_covariates,non_covid_deaths %>% transmute(death=as.factor(death),eid=eid),by="eid"),ICU_total,by="eid") %>% mutate(
  death = as.factor(ifelse(is.na(death.x)&is.na(death.y),0,1)),
  bressupdays=ifelse(is.na(bressupdays),0,bressupdays),
  aressupdays=ifelse(is.na(aressupdays),0,aressupdays),
  bcardsupdays=ifelse(is.na(bcardsupdays),0,bcardsupdays),
  acardsupdays=ifelse(is.na(acardsupdays),0,acardsupdays),
  rensupdays=ifelse(is.na(rensupdays),0,rensupdays),
  neurosupdays=ifelse(is.na(neurosupdays),0,neurosupdays),
  gisupdays=ifelse(is.na(gisupdays),0,gisupdays),
  dermsupdays=ifelse(is.na(dermsupdays),0,dermsupdays),
  orgsupmax=ifelse(is.na(orgsupmax),0,orgsupmax),
  icu_stay=fct_explicit_na(icu_stay, "0"),
  bsc_res_sup=fct_explicit_na(bsc_res_sup, "0"),
  adv_res_sup=fct_explicit_na(adv_res_sup, "0"),
  other_organ_support=fct_explicit_na(other_organ_support, "0"),
  case=fct_explicit_na(as.factor(case), "0")) %>% select(-death.x,-death.y)

ICU_ukbiobank <- left_join(ICU_ukbiobank, dff2 %>% transmute(eid=eid,case=results),by="eid") %>% 
  mutate(case = as.factor(ifelse(case.x==0 & is.na(case.y),0,1))) %>% select(-case.x,-case.y)

ICU_ukbiobank_0 <-ICU_ukbiobank %>% filter(eid %in% deathdeets$eid) %>% mutate(covid_death=as.factor(1))
ICU_ukbiobank_1 <-ICU_ukbiobank %>% filter(!eid %in% deathdeets$eid) %>% mutate(covid_death=as.factor(0))
ICU_ukbiobank<-rbind(ICU_ukbiobank_1,ICU_ukbiobank_0)


##########

#cov_positive<- icu_death_ukbb%>% 
#  mutate(results=ifelse(results==1|case==1,1,0)) %>%
#  filter(results==1) %>% 
#  select(-c(case,eid)) %>%
#  mutate(results = ifelse(is.na(results),0,results)) %>% na.omit()

#### see the regressions with the Leukocytes and results they yield??

mm1<- lm(Leukocyte_telomere_length ~ leukocyte_count,ICU_ukbiobank)
sjPlot::tab_model(m1)

#leuko, see later:

# ICU regressions:
icu_general<-glm(icu_stay~Age+Sex+Etshnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% select(-case,-death))
summ_icu<-summary(icu_general)
sjPlot::tab_model(icu_general,df.method = "wald")

cov_iff_icu<-glm(case~Age+Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% filter(icu_stay==1)%>%rename(Covid_positive=case))
summ_icu_cov<-summary(cov_iff_icu)

sjPlot::tab_model(cov_iff_icu,df.method = "wald")

summ_icu<-data.frame(coef(summ_icu_cov)) %>% mutate(names=rownames(coef(summ_icu_cov)))

#################################################Death regressions:
death_general<-glm(Death_biobank~leukocyte_count+Leukocyte_telomere_length+Age+Sex+Ethnicity+Diabetes+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=binomial,
                   data=ICU_ukbiobank%>% select(-case)%>%rename(Death_biobank=death))

summ_death<-summary(death_general)

sjPlot::tab_model(death_general,df.method = "wald",p.style = "stars",terms = c("Leukocyte_telomere_length","leukocyte_count","Age","SexMale","Smoking_statusYes","DiabetesYes"))


cov_iff_death<-glm(Covid_positive~leukocyte_count+Leukocyte_telomere_length+Age+Sex+Ethnicity+Diabetes+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,
                   family=quasibinomial,
                   data=ICU_ukbiobank%>% filter(death==1)%>%rename(Covid_positive=case))

summ_death<-summary(cov_iff_death)
summ_death<- data.frame(coef(summ_death)) %>% mutate(names=rownames(coef(summ_death)))

sjPlot::tab_model(cov_iff_death,df.method = "wald",p.style = "stars",terms = c("leukocyte_count", "Age","SexMale","Vascular_diseaseHBP","Smoking_statusYes","DiabetesYes","EthnicityAsian/Asian British"))

sjPlot::tab_model(
   cov_iff_death,death_general,
  df.method = "wald",
  p.style = "stars",
  collapse.ci = TRUE,
  CSS = list(
    css.centeralign = 'text-align: left;', 
    css.firsttablecol = 'font-weight: bold;', 
    css.summary = 'color: grey;',
    css.firsttablerow = 'font-weight: bold; font-style:normal;'),
  dv.labels = c("Covid positive", "Death in the biobank"),
  terms = c("leukocyte_count","Leukocyte_telomere_length","Age","SexMale","Vascular_diseaseHBP","Smoking_statusYes","DiabetesYes","EthnicityAsian/Asian British"))


#Maybe attempt a lasso on this? and end.

ggplot(summ_death,aes(x=reorder(names,-Estimate),y=exp(Estimate)))+
  geom_point() +
  ggtitle("Death relative Ukbiobank and Covid cases")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), width=0)+
  xlab("Variables and levels") +
  ylab("OR for relative death vs death in the biobank")

ggplot(summ_icu,aes(x=reorder(names,-Estimate),y=exp(Estimate)))+
  xlab("Variables and levels") +
  ylab("OR for relative ICU admissions vs historical admissions in ICU in the biobank")+
  ggtitle("Icu Admissions")+
  geom_point() + 
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), width=0)
  
saveRDS(ICU_ukbiobank,"/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/icu_biobank.rds")
