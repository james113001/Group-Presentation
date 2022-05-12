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

dff3<- readRDS("DATASET/data_ukbiobank.rds") %>%
  mutate(myocardial_date=as.Date(myocardial_date,"%Y-%m-%d"),
         cancer_date = as.Date(cancer_date,"%Y-%m-%d"),
         stdate_dementia=as.Date(stdate_dementia,"%Y-%m-%d"), Location = as.factor(Location)) %>% select(-c(myocardial_date,stdate_dementia,cancer_date,eid))


###########
### LTL

#break it down so it's easier to run:

#See "covid_plot_univar.sh" to see how to submit the file with the chunks data

#break into chunks or silence if necessary:
dff4<-dff3 %>% select(-c(Z_adj_TS))
ids=as.character(cut(1:length(dff4)-1, breaks = nchunks, labels = 1:nchunks))
pvalues=data.frame(dff4[,ids == ichunk] , dff3 %>% select(c(Z_adj_TS, Age_recruited, Sex))) 


#divide factors and continuous:
Xvar_cts<-data.frame(pvalues[,!sapply(pvalues,is.factor)] %>% select(-Z_adj_TS))
Xvar_factors<- data.frame(pvalues[,sapply(pvalues,is.factor)])

#Generalized likelihood ratio test: is variable significant?
lrt__a<- function(Xvar) {
  mod<-summary(lm(Z_adj_TS ~ Xvar + Age_recruited + Sex , data=pvalues))
  LRT_ltl<- pchisq(mod$null-mod$deviance, mod$df.null-mod$df.residual,lower.tail = F)
  return(LRT_ltl)
}

#Check each variable one by one to see which significant:

Betas_pvals<-function(X){
  model1 = summary(lm(Z_adj_TS ~ X + Age_recruited + Sex, data=pvalues))$coefficients
  Beta = model1[2:nrow(model1),"Estimate"]
  Pvalues = model1[2:nrow(model1),"Pr(>|t|)"]
  varbl_names = rownames(model1)[2:nrow(model1)]
  df<-data.frame("Variables"=varbl_names,"beta"=Beta,"pvals"=Pvalues)
  return(df)
}

cts<-apply(Xvar_cts,2,FUN=Betas_pvals)
if (length(Xvar_factors)!=0) {factors<-apply(Xvar_factors,2,FUN=Betas_pvals)}

ltl_results<-rbind(if (length(Xvar_factors)!=0) {do.call(rbind, Map(data.frame, factors))},do.call(rbind, Map(data.frame, cts)))

ifelse(dir.exists("MultiResults_ltl"),"",dir.create("MultiResults_ltl"))
saveRDS(data.frame(ltl_results), paste0("../ANALYSIS_STATISTICS/univariate_ltl_ukbb/MultiResults_ltl/multivariate_ltl", ichunk, ".rds"))

