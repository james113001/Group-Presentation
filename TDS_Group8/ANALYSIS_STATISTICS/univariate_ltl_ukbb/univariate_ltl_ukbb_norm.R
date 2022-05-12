#
rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

######

args <- commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

## DATA LOADING

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")

library(ggplot2)
library(dplyr)

UKbiobank<- readRDS("DATASET/data_ukbiobank.rds") %>%
  mutate(myocardial_date=as.Date(myocardial_date,"%Y-%m-%d"),
         cancer_date = as.Date(cancer_date,"%Y-%m-%d"),
         stdate_dementia=as.Date(stdate_dementia,"%Y-%m-%d"), Location = as.factor(Location)) %>% select(-c(myocardial_date,stdate_dementia,cancer_date,eid))

#Scale continuous variables - what does bin mean?
UKbiobank$IMD_town = scale(UKbiobank$IMD_town)
UKbiobank$Age_recruited = scale(UKbiobank$Age_recruited)
UKbiobank$num_cancers = scale(UKbiobank$num_cancers)
UKbiobank$SleepDur = scale(UKbiobank$SleepDur)
UKbiobank$breadcons = scale(UKbiobank$breadcons)
UKbiobank$tea_intake = scale(UKbiobank$tea_intake)
UKbiobank$Age_HBP = scale(UKbiobank$Age_HBP)
UKbiobank$BMI = scale(UKbiobank$BMI)
UKbiobank$NO2_poll = scale(UKbiobank$NO2_poll)
UKbiobank$Height = scale(UKbiobank$Height)
UKbiobank$DiaBP = scale(UKbiobank$Dia_BP)
UKbiobank$numb_household = scale(UKbiobank$numb_household)
UKbiobank$IMD_ = scale(UKbiobank$IMD)
UKbiobank$heart_attack_bin = scale(UKbiobank$heart_attack_bin)
UKbiobank$medic_comorbidities_bin = scale(UKbiobank$medic_comorbidities_bin)
UKbiobank$Age_diabetes_bin = scale(UKbiobank$Age_diabetes_bin)
UKbiobank$length_current_address = scale(UKbiobank$length_current_address)
UKbiobank$num_days_walk_10mins = scale(UKbiobank$num_days_walk_10mins)
UKbiobank$cooked_veggie_intake = scale(UKbiobank$cooked_veggie_intake)
UKbiobank$fresh_fruit_intake = scale(UKbiobank$fresh_fruit_intake)
UKbiobank$vit_suplements = scale(UKbiobank$vit_suplements)
UKbiobank$disability_allowance = scale(UKbiobank$disability_allowance)
UKbiobank$NoEDWS = scale(UKbiobank$NoEDWS)
UKbiobank$common_disease_diag = scale(UKbiobank$common_disease_diag)
UKbiobank$mon_count_dev = scale(UKbiobank$mon_count_dev)


###########
### LTL

#break it down so it's easier to run:

#See "covid_plot_univar.sh" to see how to submit the file with the chunks data

#break into chunks or silence if necessary:
dff5<-UKbiobank %>% select(-Z_adj_TS)
ids=as.character(cut(1:length(dff5)-1, breaks = nchunks, labels = 1:nchunks))
pvalues=data.frame(dff5[,ids == ichunk] , "Z_adj_TS"= UKbiobank %>% 
                     select(Z_adj_TS) %>% 
                     mutate(Z_adj_TS = as.numeric(Z_adj_TS)))


#divide factors and continuous:
Xvar_cts<-data.frame(pvalues[,!sapply(pvalues,is.factor)] %>% select(-Z_adj_TS))
Xvar_factors<- data.frame(pvalues[,sapply(pvalues,is.factor)])

#Generalized likelihood ratio test: is variable significant?
lrt__a<- function(Xvar) {
  mod<-summary(lm(Z_adj_TS ~ Xvar , UKbiobank=pvalues))
  LRT_ltl<- pchisq(mod$null-mod$deviance, mod$df.null-mod$df.residual,lower.tail = F)
  return(LRT_ltl)
}

#Check each variable one by one to see which significant:

Betas_pvals<-function(X){
  model1 = summary(lm(Z_adj_TS ~ X, data=pvalues))$coefficients
  Beta = model1[2:nrow(model1),"Estimate"]
  Pvalues = model1[2:nrow(model1),"Pr(>|t|)"]
  varbl_names = rownames(model1)[2:nrow(model1)]
  df<-data.frame("Variables"=varbl_names,"beta"=Beta,"pvals"=Pvalues)
  return(df)
}

cts<-apply(Xvar_cts,2,FUN=Betas_pvals)
if (length(Xvar_factors)!=0) {factors<-apply(Xvar_factors,2,FUN=Betas_pvals)}

ltl_results<-rbind(if (length(Xvar_factors)!=0) {do.call(rbind, Map(data.frame, factors))},do.call(rbind, Map(data.frame, cts)))

ifelse(dir.exists("Results_ltl_normalised"),"",dir.create("Results_ltl_normalised"))
saveRDS(data.frame(ltl_results), paste0("../ANALYSIS_STATISTICS/univariate_ltl_ukbb/Results_ltl/univariate_ltl_norm", ichunk, ".rds"))

