
setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outputs")
biodata<-readRDS("official_dataset.rds")

library(ggplot2)
library(dplyr)

colnames(biodata[c(30:60)])
summary(biodata[c(30:60)],na.rm=TRUE)

colnames(biodata1)

hist(biodata$salt_in_food)

biodata1<-biodata %>% 
    select(-age_st_smoking,-father_age,-Age_HBP,-Age_diabetes) %>%  #remove variables with > 75% NA 
    mutate(tobacco_exp_outside_hshld_up=ifelse(is.na(tobacco_exp_outside_hshld), 0, tobacco_exp_outside_hshld)) #%>% 
    #mutate(to_filter = ifelse(is.na(tobacco_exp_outside_hshld), 0,tobacco_exp_outside_hshld)) %>% 

)

summary(biodata1["tobacco_exp_outside_hshld_up"],na.rm=TRUE)
