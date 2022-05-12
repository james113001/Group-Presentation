
library(ggplot2)
library(dplyr)

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outputs")
biobank<-readRDS("official_dataset.rds")

colnames(biobank[c(60:90)])
summary(biobank[c(60:90)],na.rm=TRUE)


# Data Distributions 
par(mfrow=c(4,2))
hist(biobank$Age_diabetes, main='Age diabetes', col="lightblue");
hist(biobank$forces_vital_capacity, main = "forces_vital_capacity", col = "lightblue");
hist(biobank$MenopauseAge, main = "MenopauseAge", col = "lightblue");
hist(biobank$age_heart_attack, main = "age_heart_attack",col = "lightblue");
hist(biobank$medic_comorbidities, main = "medic_comorbidities", col = "lightblue");
hist(biobank$vit_suplements, main = "disability_allowance", col = "lightblue")
par(mfrow=c(4,2))
hist(biobank$neuroticism_score, main = "neuroticism_score", col = "lightblue")
hist(biobank$BMI, main = "BMI", col = "lightblue")
hist(biobank$NoEDWS, main = "NoEDWS", col = "lightblue")
hist(biobank$common_disease_diag, main = "common_disease_diag", col = "lightblue")
hist(biobank$Age_diabetes, main = "Age_diabetes", col = "lightblue")
hist(biobank$age_heart_attack, main = "age_heart_attack", col = "lightblue")


# Remove outliers for forces vital capacity?
summary(biobank$forces_vital_capacity)
boxplot(biobank$forces_vital_capacity)


biobank$heart_attack_YN <- ifelse(!is.na (biobank$age_heart_attack), 1, 0) 
biobank$medic_comorbidities[is.na(biobank$medic_comorbidities)] <- 0
View(biobank$medic_comorbidities)

View(biobank)

## Izzy edits: Creating a heart attack Yes/No column and recoding NAs in medic comorbidities, diabetes age, heart attack age and menopause age to 0 
data$heart_attack_YN <- ifelse(!is.na (data$age_heart_attack), 1, 0) 
data$medic_comorbidities[is.na(data$medic_comorbidities)] <- 0
data$Age_diabetes[is.na(data$Age_diabetes)] <- 0
data$age_heart_attack[is.na(data$age_heart_attack)] <- 0
data$MenopauseAge[is.na(data$MenopauseAge)] <- 0

# Very High missingness
# Keep: Age of diabetes, age heart attack - keep because NA could indicate no diabetes/heartattack


# Delete: num of still births, num of miscarriages, COPD diag, weed freq, Bipolar depression, Addicted drugs, EBV, CMV, HBV  - Remove because of very high NA?

#  High missing-ness
# Keep: Age of menopause

#  Mid
# Medic comorbidities - NA --> 0



hist(biobank$illnesses_father, main = "illnesses_father", col = "lightblue")
hist(biobank$illnesses_mother, main = "illnesses_mother", col = "lightblue")
hist(biobank$illnesses_siblings, main = "illnesses_siblings", col = "lightblue")
hist(biobank$smoking_status, main = "smoking_status", col = "lightblue")
hist(biobank$alcohol_drink_status, main = "alcohol_drink_status", col = "lightblue")
hist(biobank$Cannabis, main = "Cannabis", col = "lightblue")
hist(biobank$feeling_depression, main = "feeling_depression", col = "lightblue")
hist(biobank$IBS, main = "IBS", col = "lightblue")
hist(biobank$feeling_depression, main = "feeling_depression", col = "lightblue")




#From NEL: Decide imputation/removal/ignoring - MNAR/MCAR/other
Data<-data[1:1000,]
na_perc<-round(100*sapply(data, function(x) sum(length(which(is.na(x)))))/length(data$Sex),3)
rownames(data.frame(na_perc))

#Good variables
data.frame(na_perc[na_perc<25])

#Middling variables
na_perc[na_perc<50 &na_perc>25]

#poor Quality variables
mid<-data.frame(na_perc[na_perc>50 & na_perc<=70])
mid
low<-data.frame(na_perc[na_perc>70])
low
filter_missings <- Data %>% select(-c(rownames(low),rownames(mid)))




########################
