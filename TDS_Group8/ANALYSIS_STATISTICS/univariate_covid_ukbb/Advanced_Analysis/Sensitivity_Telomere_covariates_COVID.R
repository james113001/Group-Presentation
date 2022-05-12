library(dplyr)
library(glmnet)
library(focus)
library(tidyverse)
library(forcats)
library(withr)

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis")

#Get the dataset with the extra varables, but more constrained:

ukbiobank_dataset<-readRDS("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data/DATASET/trees_dataset.rds")

#Add:: Ethnicity and cycstatin, what happened with Alcohol_drinker_status

Covid_covariates<-ukbiobank_dataset%>% 
  select("eid","Age","Ethnicity","Sex","Diabetes","High_blood_pressure","Vascular_disease","Average_household_income","leukocyte_count","Myocardial_disease","Smoking_status",
         "Alcohol_drinker_status","Leukocyte_telomere_length","Walk_pace","Fresh_fruit_intake","Oily_fish_intake","Mother_smoking_birth",
         "Processed_meat_consumption","Disability_allowance","Index_of_multiple_deprivation","BMI")%>%
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

colnames(icu_death_ukbb)

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

setwd("/rds/general/project/hda_21-22/live/TDS/General/Data")
deathcause<-fread("death_cause.txt")
deathdates<-fread("death.txt")

deathdeets<-deathcause %>%
  filter(cause_icd10=="U071")%>%
  mutate(eid=as.character(eid)) %>% distinct() %>% select(eid)


ICU_ukbiobank_0 <-ICU_ukbiobank %>% filter(eid %in% deathdeets$eid) %>% mutate(covid_death=as.factor(1))
ICU_ukbiobank_1 <-ICU_ukbiobank %>% filter(!eid %in% deathdeets$eid) %>% mutate(covid_death=as.factor(0))
ICU_ukbiobank<-rbind(ICU_ukbiobank_1,ICU_ukbiobank_0)

##########

#### see the regressions with the Leukocytes and results they yield??

#library(heatmap)
summary(ICU_ukbiobank$covid_death)

# Death regressions:

y <- ICU_ukbiobank%>% filter(death==1) %>% select(covid_death) %>% mutate(covid_death=as.numeric(covid_death)-1)
x <- ICU_ukbiobank%>% filter(death==1) %>% select(-covid_death,-case,-eid) 

x_train <- model.matrix( ~ .-1, x)
glm<-cv.glmnet(x = x_train, y = as.matrix(y), penalty.factor=c(rep(1,ncol(x_train))), family = "binomial")

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_glmnet_death.pdf"  ))
plot(glm)
dev.off()

beta_lasso=coef(glm, s="lambda.min")[2:(ncol(x_train)+1),]
selected_lasso=names(beta_lasso)[which(beta_lasso!=0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

out = VariableSelection(x = x_train,seed=1,y = as.matrix(y),family = "binomial")

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_calib_death.pdf"))
CalibrationPlot(out)
dev.off()

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)

selection_plot<-data.frame(names=names(selprop),prop=selprop)%>%mutate(good=as.factor(ifelse(selprop>hat_params[2],"1","0"))) %>% arrange(desc(prop))

a <- ifelse(selection_plot$good == "0", "#999999", "#3c7be8")
pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_heatmap_death.pdf"  ))
ggplot(selection_plot,aes(x=fct_reorder(names,desc(prop)),y=prop,col=good,fill=good))+
  geom_segment(aes(yend=0,xend=names),size=2,lineend = "round")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1,colour=a),
        text = element_text(size=13),
        strip.placement = "outside",
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))



plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}

dev.off()

# ICU
removal_icu<-c("bressupdays", "aressupdays", "bcardsupdays", "acardsupdays", "rensupdays", "neurosupdays", "gisupdays", "dermsupdays","orgsupmax","bsc_res_sup","adv_res_sup","other_organ_support","death")


y<- ICU_ukbiobank%>% filter(icu_stay==1) %>% select(case) %>% mutate(case=as.numeric(case)-1)
x<- ICU_ukbiobank%>% filter(icu_stay==1) %>% select(-c(case,eid))  %>% select(-removal_icu)

x_train <- model.matrix( ~ .-1, x)
glm<-cv.glmnet(x = x_train, y = as.matrix(y), penalty.factor=c(rep(1,ncol(x_train))), family = "binomial")

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_glmnet_icu.pdf"  ))
plot(glm)
dev.off()

beta_lasso=coef(glm, s="lambda.min")[2:(ncol(x_train)+1),]
selected_lasso=names(beta_lasso)[which(beta_lasso!=0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

out = VariableSelection(x = x_train,seed=1,y = as.matrix(y),family = "binomial")

#pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_calib_icu.pdf"))
#CalibrationPlot(out)
#dev.off()

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

selection_plot<-data.frame(names=names(selprop),prop=selprop)%>%mutate(good=as.factor(ifelse(selprop>hat_params[2],"1","0"))) %>% arrange(desc(prop))

a <- ifelse(selection_plot$good == "0", "#999999", "#3c7be8")
png(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/Sensitivity_plots/plot_heatmap_death.png"))
ggplot(selection_plot,aes(x=fct_reorder(names,desc(prop)),y=prop,col=good,fill=good))+
  geom_segment(aes(yend=0,xend=names),size=2,lineend = "round")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1,colour=a),
        text = element_text(size=13),
        strip.placement = "outside",
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))

dev.off()

plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}




#Maybe attempt a lasso on this? and end.
