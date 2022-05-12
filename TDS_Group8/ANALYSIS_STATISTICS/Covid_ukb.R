


#########
#      OBSOLETE FILE - SEE "agile_apply_coding" for the better version!!
########

# Combine the datasets and the variables 
rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")

library(ggplot2)
library(dplyr)

dff1<-readRDS("DATASET/deaths_icu_covid.rds")%>% 
  mutate(eid = as.character(eid),
         Loc=as.factor(Loc)) %>% 
  select(-liversupdays) 

dff2<-readRDS("DATASET/test_covid.rds") %>% mutate(results=1, eid = as.character(eid))
dff3<- readRDS("DATASET/data_ukbiobank.rds") %>%
  mutate(myocardial_date=as.Date(myocardial_date,"%Y-%m-%d"),
         cancer_date = as.Date(cancer_date,"%Y-%m-%d"),
         stdate_dementia=as.Date(stdate_dementia,"%Y-%m-%d")) %>% select(-c(myocardial_date,stdate_dementia,cancer_date))

test_ukbb<-left_join(dff3,dff2) %>% mutate(results = ifelse(is.na(results),0,results))

icu_death_ukbb<-left_join(left_join(dff3,dff1) %>% mutate(death = as.factor(ifelse(is.na(death),0,death)),
                                                          ccadmitype=ifelse(is.na(ccadmitype),0,ccadmitype),
                                                          ccadmisorc=ifelse(is.na(ccadmisorc),0,ccadmisorc),
                                                          bressupdays=ifelse(is.na(bressupdays),0,bressupdays),
                                                          aressupdays=ifelse(is.na(aressupdays),0,aressupdays),
                                                          bcardsupdays=ifelse(is.na(bcardsupdays),0,bcardsupdays),
                                                          acardsupdays=ifelse(is.na(acardsupdays),0,acardsupdays),
                                                          rensupdays=ifelse(is.na(rensupdays),0,rensupdays),
                                                          neurosupdays=ifelse(is.na(neurosupdays),0,neurosupdays),
                                                          icu_stay = as.factor(ifelse(is.na(ccstartdate),0,1))),dff2) %>% 
                                                          select(-result) %>%
                                                          mutate(results = ifelse(is.na(results),0,results))


################################
### univariate analysis of the biobank, against covid ICU data work with chunks:

# count number of levels in the variables:
max_levels<-1

exclusion_vars_icu<-colnames(icu_death_ukbb)[c(53,57:62,73:78,82)]
data<-icu_death_ukbb %>% select(-exclusion_vars_icu) %>% filter(results==1)

for (i in 1:length(data)){
  max_levels<-  ifelse(length(levels(data[,i]))<max_levels,max_levels,
                       length(levels(data[,i])))
}

Beta = Pvalues = Names= matrix(NA, nrow = max_levels,ncol = (length(data)-1))

## split this into different nodes and run better
ids=as.character(cut(1:length(data), breaks = nchunks, labels = 1:nchunks))

pvalues=data.frame(cbind(data[,ids == ichunk ],icu_stay=data$icu_stay))
GLRT_icu<-matrix(NA, nrow = ncol(pvalues),ncol = 1)

print(colnames(pvalues)[1:length(pvalues)])
summary(pvalues)

for (j in  1:(length(pvalues)-1)) {
  print(paste0("errors_in_glrt",j))
  mod<-summary(glm(icu_stay ~pvalues[,j] , data=pvalues,family=binomial))
  GLRT_icu[j,1]<- pchisq(mod$null-mod$deviance, mod$df.null-mod$df.residual,lower.tail = F)
}


for (j in  1:(length(pvalues)-1)) {
  # print(j) debugging tool
  if(class(pvalues[,j])=="factor" & colnames(pvalues[j])!="icu_stay"){
    len<-length(levels(pvalues[,j]))
    model1 = glm(icu_stay ~ pvalues[,j], data=pvalues,family=binomial)
    Beta[, j] = c(coefficients(model1)[2:c(len)],rep(NA,max_levels-len+1))
    Names[, j] = c(paste0(colnames(pvalues[j])," ",row.names(summary(model1)$coefficients)[2:c(len)]),rep(NA,max_levels-len+1))
    Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                     rep(NA,max_levels+1-c(length(summary(model1)$coefficients[,1]))))
  }
  if(class(pvalues[,j]) != "factor" & colnames(pvalues[j])!="icu_stay"){
    model1 = glm(icu_stay ~ pvalues[,j], data=pvalues,family=binomial)
    Beta[, j] = c(coefficients(model1)[2],rep(NA,max_levels-1))
    Names[, j] = c(colnames(pvalues[j]),rep(NA,max_levels-1))
    Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                     rep(NA,max_levels+1-c(length(summary(model1)$coefficients[,1]))))
  }
  print(paste0("no_error_icu:",colnames(pvalues)[j]))}

plot_betas<-data.frame("var_names"=c(Names),"p_val"=as.numeric(Pvalues),
                       "beta"=c(Beta)) %>% filter(is.na(var_names)==FALSE)

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb")

ifelse(dir.exists("Results"),"",dir.create("Results"))
saveRDS(plot_betas, paste0("Results/univ_pvalues_chnk_", ichunk, ".rds"))

################################
### univariate analysis of the biobank, against covid-death data work with chunks:
# count number of levels in the variables:

max_levels<-1
death_excl_vars<-colnames(icu_death_ukbb)[c(53,57:62,74:78,80)]
data<-icu_death_ukbb %>% select(-death_excl_vars) %>% filter(results==1)

for (i in 1:length(data)){
  max_levels<-  ifelse(length(levels(data[,i]))<max_levels,max_levels,length(levels(data[,i])))
}

Beta = Pvalues = Names= matrix(NA, nrow = max_levels,ncol = (length(data)-1))

## split this into different nodes and run better

ids=as.character(cut(1:length(data), breaks = nchunks, labels = 1:nchunks))

pvalues=data.frame(cbind(data[,ids == ichunk ],death=data$death))

GLRT_death<-matrix(NA, nrow = ncol(pvalues),ncol = 1)
colnames(pvalues)

for (j in  1:(length(pvalues)-1)) {
  mod<-summary(glm(death ~ pvalues[,j], data=pvalues,family=binomial))
  GLRT_death[j,1]<-  pchisq(mod$null-mod$deviance, mod$df.null-mod$df.residual,lower.tail = F)
}

for (j in  1:(length(pvalues)-1)) {
  # print(j) debugging tool
  if(class(pvalues[,j])=="factor" & colnames(pvalues[j])!="death"){
    
    len<-length(levels(pvalues[,j]))
    
    model1 = glm(death ~ pvalues[,j], data=pvalues,family=binomial)
    Beta[, j] = c(coefficients(model1)[2:c(len)],rep(NA,max_levels-len+1))
    Names[, j] = c(paste0(colnames(pvalues[j])," ",row.names(summary(model1)$coefficients)[2:c(len)]),rep(NA,max_levels-len+1))
    Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                     rep(NA,max_levels+1-c(length(summary(model1)$coefficients[,1]))))
  }
  
  if(class(pvalues[,j]) !="factor" & colnames(pvalues[j])!="death"){
    model1 = glm(death ~ pvalues[,j], data=pvalues,family=binomial)
    Beta[, j] = c(coefficients(model1)[2],rep(NA,max_levels-1))
    Names[, j] = c(colnames(pvalues[j]),rep(NA,max_levels-1))
    Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                     rep(NA,max_levels+1-c(length(summary(model1)$coefficients[,1]))))
  }
  print(paste0("no_error_death:",colnames(pvalues)[j]))
  }

plot_betas<-data.frame("var_names"=c(Names),"p_val"=as.numeric(Pvalues),
                       "beta"=c(Beta)) %>% filter(is.na(var_names)==FALSE)

ifelse(dir.exists("Results_death"),"",dir.create("Results_death"))
saveRDS(plot_betas, paste0("./Results_death/univ_death_pval_chnk_", ichunk, ".rds"))

gltr<-data.frame(deaths_glrt=GLRT_death,icu_glrt=GLRT_icu)

ifelse(dir.exists("Results_gltr"),"",dir.create("Results_gltr"))
saveRDS(gltr, paste0("./Results_gltr/gltr", ichunk, ".rds"))

