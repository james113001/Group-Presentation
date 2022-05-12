rm(list=ls())
library(glmnet)
library(skimr)
library(broom)
library(ggplot2)
library(dplyr)
library(janitor)
library(viridis)
library(focus)
library(tidyverse)

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")
#dataLTL<-readRDS("DATASET/deaths_icu_covid.rds")
ukbiobank_subgrouped<-readRDS("DATASET/ukbiobank_subgrouped.rds")
dff<- readRDS("DATASET/data_ukbiobank_updated.rds")

LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}
LoadPackages(c("focus","igraph","glmnet", "pheatmap"))


#Lasso Regression on biobank data for LTL 

dff<- data.frame(dff) %>% select(colnames(dff) %in% c(names_variable[,1]))

ltl <- dff %>% 
  dplyr::select(Leukocyte_telomere_length) 
ltlpred<-dff %>% 
  dplyr::select(-Leukocyte_telomere_length, -eid)
#ltlpred<-data.matrix(ltlpred)
ltlpred<-model.matrix( ~ ., data.frame(ltlpred))

names_variable<-ukbiobank_subgrouped%>% na.omit() %>% select(variables)
empty_list<-matrix(NA,nrow=2*length(colnames(ltlpred))+2,ncol=2)
z=1

for (i in 1:length(names_variable[,1])){
  hold_this_beer<-data.frame(var_names=colnames(ltlpred)) %>% mutate(var_names=as.character(var_names)) %>% filter(str_detect(var_names, names_variable[i,1]))
  hold_this_beer
  for (j in 1:length(hold_this_beer[,'var_names'])){
    empty_list[z,1]<-names_variable[i,1]
    empty_list[z,2]<-if (length(hold_this_beer[j,1])>0){hold_this_beer[j,1]} else (NA)
    z=z+1  
  }
  print(i)
}

subgrouped_list<-left_join(data.frame(empty_list %>% na.omit()) %>% rename(variables=X1,levels=X2),ukbiobank_subgrouped)
subgrouped_list<-subgrouped_list[-120,]

#tukbiobank_subgrouped<- data.frame(t(ukbiobank_subgrouped))%>% na.omit() %>% row_to_names(1)
#subgroups <- tukbiobank_subgrouped %>% select(colnames(ltlpred))

ltlfem <- dff %>% 
  filter(Sex == 'Female') %>%
  dplyr::select(Leukocyte_telomere_length) 
ltlpredfem <- dff %>% 
  filter(Sex == 'Female') %>%
  dplyr::select(-Leukocyte_telomere_length, -eid)
ltlpredfem<-model.matrix( ~ ., data.frame(ltlpredfem))



ltlcv_model <- cv.glmnet(ltlpred, ltl$Leukocyte_telomere_length, alpha = 1)
ltlbest_lambda <- ltlcv_model$lambda.min
ltlbest_lambda
#2.87E-4

ltlbest_model <- glmnet(ltlpred, ltl$Leukocyte_telomere_length, alpha = 1, lambda = ltlbest_lambda)
ltltidy_best_model <- tidy(ltlbest_model)
coef(ltlbest_model)
write.csv(ltltidy_best_model, "tidyltl_best_model.csv")

#Stability selection
out=VariableSelection(xdata=ltlpred, ydata=ltl$Leukocyte_telomere_length, verbose=FALSE, 
                      penalty.factor=c(rep(1,ncol(ltlpred))),
                      family="gaussian")

png(file = "/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data/cali.png")
CalibrationPlot(out)
dev.off()

###############################females
ltlfemcv_model <- cv.glmnet(ltlpredfem, ltlfem$Leukocyte_telomere_length, alpha = 1)
ltlfembest_lambda <- ltlfemcv_model$lambda.min
ltlfembest_lambda
#9.21E-4

ltlfembest_model <- glmnet(ltlpredfem, ltlfem$Leukocyte_telomere_length, alpha = 1, lambda = ltlfembest_lambda)
ltlfemtidy_best_model <- tidy(ltlfembest_model)
coef(ltlfembest_model)
write.csv(ltlfemtidy_best_model, "tidyfemltl_best_model.csv")

#Stability selection
outfem=VariableSelection(xdata=ltlpredfem, ydata=ltlfem$Leukocyte_telomere_length, verbose=FALSE, 
                      penalty.factor=c(rep(1,ncol(ltlpredfem))),
                      family="gaussian")
png(file = "/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data/fem_cali.png")
CalibrationPlot(outfem)
dev.off()

# Calibrated selection proportions 
selprop=SelectionProportions(out)
selpropfem=SelectionProportions(outfem)


# Calibrated parameters
hat_params=Argmax(out)
hat_paramsfem=Argmax(outfem)


##########SELECTION PROPORTION GRAPH###########################################################
###############################################################################################

hmm<-data.frame(Feature=names(selprop), Sel_Prop=selprop, row.names=NULL)
hmm<- left_join(hmm %>% filter(Feature!="(Intercept)"),subgrouped_list %>%mutate(Feature=levels)%>% select(type,Feature))
hmm$good<-as.factor(ifelse(hmm$Sel_Prop>hat_params[2],"1","0"))

hmmfem<-data.frame(Feature=names(selpropfem), Sel_Prop=selpropfem, row.names=NULL)
hmmfem<- left_join(hmmfem %>% filter(Feature!="(Intercept)"),subgrouped_list %>%mutate(Feature=levels)%>% select(type,Feature))
hmmfem$good<-as.factor(ifelse(hmmfem$Sel_Prop>hat_paramsfem[2],"1","0"))
hmmfem<- hmmfem[-1,] ##remove Sex



ggplot(hmm, aes(x=fct_reorder(Feature,desc(Sel_Prop)), y=Sel_Prop, color = good, fill =good))+
  facet_grid(~hmm$type, scales = "free_x", space = "free_x")+ 
  geom_segment(aes(yend=0,xend=Feature),size=2,lineend = "round")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion")+
  xlab("")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
  text = element_text(size=10),
  strip.placement = "outside",
  legend.position="none",
  strip.background = element_rect(fill = "white"),
  plot.title = element_text(size=14, face="bold.italic"))

ggplot(hmmfem, aes(x=fct_reorder(Feature,desc(Sel_Prop)), y=Sel_Prop, color = good, fill =good))+
  facet_grid(~hmmfem$type, scales = "free_x", space = "free_x")+ 
  geom_segment(aes(yend=0,xend=Feature),size=2,lineend = "round")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion for Females Only")+
  xlab("")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
                    text = element_text(size=10),
                    strip.placement = "outside",
                    legend.position="none",
                    strip.background = element_rect(fill = "white"),
                    plot.title = element_text(size=14, face="bold.italic"))