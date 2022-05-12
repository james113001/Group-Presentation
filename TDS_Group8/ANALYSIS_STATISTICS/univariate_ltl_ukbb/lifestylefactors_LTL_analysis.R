getwd()

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data/DATASET")
library(ggplot2)
library(dplyr)
library(tidyverse)
suppressPackageStartupMessages(library(tidyr))
library(devtools)
install.packages('glmnet', version = "2.0_5")

install.packages('glmnet', version = "4.1-3")
library(glmnet)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) 
  dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) 
  file.create(M)
cat("\nCXX14FLAGS=-O3 -Wno-unused-variable -Wno-unused-function",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -Wno-unused-variable -Wno-unused-function",
    file = M, sep = "\n", append = TRUE)


#####Univariate analysis to find variables driving LTL - use Nels code




#####Multivariate analysis to find variables driving LTL
UKbiobank <- read_rds("data_ukbiobank.rds") 

#Intgers into num
UKbiobank$length_current_address <- as.numeric(UKbiobank$length_current_address)
UKbiobank$num_days_walk_10mins <- as.numeric(UKbiobank$num_days_walk_10mins)
UKbiobank$cooked_veggie_intake <- as.numeric(UKbiobank$cooked_veggie_intake)
UKbiobank$fresh_fruit_intake <- as.numeric(UKbiobank$fresh_fruit_intake)
UKbiobank$vit_suplements <- as.numeric(UKbiobank$vit_suplements)
UKbiobank$disability_allowance <- as.numeric(UKbiobank$disability_allowance)
UKbiobank$NoEDWS <- as.numeric(UKbiobank$NoEDWS)
UKbiobank$common_disease_diag <- as.numeric(UKbiobank$common_disease_diag)
UKbiobank$mon_count_dev <- as.numeric(UKbiobank$mon_count_dev)


#Scale continuous variables - what does bin mean?
UKbiobank$IMD_town_normalized = scale(UKbiobank$IMD_town)
UKbiobank$Age_recruited_normalized = scale(UKbiobank$Age_recruited)
UKbiobank$num_cancers_normalized = scale(UKbiobank$num_cancers)
UKbiobank$SleepDur_normalized = scale(UKbiobank$SleepDur)
UKbiobank$breadcons_normalized = scale(UKbiobank$breadcons)
UKbiobank$tea_intake_normalized = scale(UKbiobank$tea_intake)
UKbiobank$Age_HBP_normalized = scale(UKbiobank$Age_HBP)
UKbiobank$BMI_normalized = scale(UKbiobank$BMI)
UKbiobank$NO2_poll_normalized = scale(UKbiobank$NO2_poll)
UKbiobank$Height_normalized = scale(UKbiobank$Height)
UKbiobank$DiaBP_normalized = scale(UKbiobank$Dia_BP)
UKbiobank$numb_household_normalized = scale(UKbiobank$numb_household)
UKbiobank$IMD_normalized = scale(UKbiobank$IMD)
UKbiobank$heart_attack_bin_normalized = scale(UKbiobank$heart_attack_bin)
UKbiobank$medic_comorbidities_bin_normalized = scale(UKbiobank$medic_comorbidities_bin)
UKbiobank$Age_diabetes_bin_normalized = scale(UKbiobank$Age_diabetes_bin)
UKbiobank$length_current_address_normalized = scale(UKbiobank$length_current_address)
UKbiobank$num_days_walk_10mins_normalized = scale(UKbiobank$num_days_walk_10mins)
UKbiobank$cooked_veggie_intake_normalized = scale(UKbiobank$cooked_veggie_intake)
UKbiobank$fresh_fruit_intake_normalized = scale(UKbiobank$fresh_fruit_intake)
UKbiobank$vit_suplements_normalized = scale(UKbiobank$vit_suplements)
UKbiobank$disability_allowance_normalized = scale(UKbiobank$disability_allowance)
UKbiobank$NoEDWS_normalized = scale(UKbiobank$NoEDWS)
UKbiobank$common_disease_diag_normalized = scale(UKbiobank$common_disease_diag)
UKbiobank$mon_count_dev_normalized = scale(UKbiobank$mon_count_dev)

#Run model- what to do with date and locatcion (Chr)
Model1 <- lm(UKbiobank$Z_adj_TS ~ UKbiobank$Sex + UKbiobank$IMD_town_normalized + 
               UKbiobank$Age_recruited_normalized + UKbiobank$num_cancers_normalized + 
               UKbiobank$SleepDur_normalized + UKbiobank$breadcons_normalized + UKbiobank$tea_intake_normalized + 
               UKbiobank$Age_HBP_normalized + UKbiobank$BMI_normalized + UKbiobank$NO2_poll_normalized + 
               UKbiobank$Height_normalized + UKbiobank$DiaBP_normalized + UKbiobank$numb_household_normalized +
               UKbiobank$IMD_normalized + UKbiobank$heart_attack_bin_normalized + UKbiobank$medic_comorbidities_bin_normalized +
               UKbiobank$Age_diabetes_bin_normalized + UKbiobank$length_current_address_normalized + UKbiobank$num_days_walk_10mins_normalized + 
               UKbiobank$cooked_veggie_intake_normalized + UKbiobank$fresh_fruit_intake_normalized + UKbiobank$vit_suplements_normalized +
               UKbiobank$disability_allowance_normalized + UKbiobank$NoEDWS_normalized + UKbiobank$common_disease_diag_normalized +
               UKbiobank$mon_count_dev_normalized + UKbiobank$accomodation_lived + UKbiobank$own_accomodation + UKbiobank$avg_household_income +
               UKbiobank$walk_pace + UKbiobank$Insomnia + UKbiobank$narcolepsy + UKbiobank$CurrSmok + UKbiobank$oily_fish_intake + UKbiobank$non_oily_fish_intake +
               UKbiobank$ProcMeat + UKbiobank$cheese_cons + UKbiobank$salt_in_food + UKbiobank$Alcohol_intake_freq + UKbiobank$mat_smoking_birth + 
               UKbiobank$dep_last_2_weeks + UKbiobank$long_standing_illness + UKbiobank$diabetes_doctor + UKbiobank$cancer_doctor + 
               UKbiobank$presc_medicine + UKbiobank$illnesses_father + UKbiobank$illnesses_mother)

summary(Model1)
library(ggplot2) 

library(broom)
Model1resultsV <- tidy(Model1)
Model1resultsV$pvalue <- as.numeric(Model1resultsV$p.value)
Model1resultsV$estimate <- as.numeric(Model1resultsV$estimate)
pval <- Model1resultsV$pvalue
View(Model1resultsV)
estimate <- Model1resultsV$estimate



#Lasso regression
LassoUKbiobank <- read_rds("data_ukbiobank.rds") 
LassoUKbiobank$length_current_address <- as.numeric(UKbiobank$length_current_address)
LassoUKbiobank$num_days_walk_10mins <- as.numeric(UKbiobank$num_days_walk_10mins)
LassoUKbiobank$cooked_veggie_intake <- as.numeric(UKbiobank$cooked_veggie_intake)
LassoUKbiobank$fresh_fruit_intake <- as.numeric(UKbiobank$fresh_fruit_intake)
LassoUKbiobank$vit_suplements <- as.numeric(UKbiobank$vit_suplements)
LassoUKbiobank$disability_allowance <- as.numeric(UKbiobank$disability_allowance)
LassoUKbiobank$NoEDWS <- as.numeric(UKbiobank$NoEDWS)
LassoUKbiobank$common_disease_diag <- as.numeric(UKbiobank$common_disease_diag)
LassoUKbiobank$mon_count_dev <- as.numeric(UKbiobank$mon_count_dev)

LassoUKbiobank$IMD_town = scale(UKbiobank$IMD_town)
LassoUKbiobank$Age_recruited = scale(UKbiobank$Age_recruited)
LassoUKbiobank$num_cancers = scale(UKbiobank$num_cancers)
LassoUKbiobank$SleepDur = scale(UKbiobank$SleepDur)
LassoUKbiobank$breadcons = scale(UKbiobank$breadcons)
LassoUKbiobank$tea_intake = scale(UKbiobank$tea_intake)
LassoUKbiobank$Age_HBP = scale(UKbiobank$Age_HBP)
LassoUKbiobank$BMI = scale(UKbiobank$BMI)
LassoUKbiobank$NO2_poll = scale(UKbiobank$NO2_poll)
LassoUKbiobank$Height = scale(UKbiobank$Height)
LassoUKbiobank$DiaBP = scale(UKbiobank$Dia_BP)
LassoUKbiobank$numb_household = scale(UKbiobank$numb_household)
LassoUKbiobank$IMD = scale(UKbiobank$IMD)
LassoUKbiobank$heart_attack_bin = scale(UKbiobank$heart_attack_bin)
LassoUKbiobank$medic_comorbidities_bin = scale(UKbiobank$medic_comorbidities_bin)
LassoUKbiobank$Age_diabetes_bin = scale(UKbiobank$Age_diabetes_bin)
LassoUKbiobank$length_current_address = scale(UKbiobank$length_current_address)
LassoUKbiobank$num_days_walk_10mins = scale(UKbiobank$num_days_walk_10mins)
LassoUKbiobank$cooked_veggie_intake = scale(UKbiobank$cooked_veggie_intake)
LassoUKbiobank$fresh_fruit_intake = scale(UKbiobank$fresh_fruit_intake)
LassoUKbiobank$vit_suplements = scale(UKbiobank$vit_suplements)
LassoUKbiobank$disability_allowance = scale(UKbiobank$disability_allowance)
LassoUKbiobank$NoEDWS = scale(UKbiobank$NoEDWS)
LassoUKbiobank$common_disease = scale(UKbiobank$common_disease_diag)
LassoUKbiobank$mon_count_dev = scale(UKbiobank$mon_count_dev)

Z_adj_TS <- which(colnames(LassoUKbiobank) == "Z_adj_TS")
Y <- as.numeric(LassoUKbiobank[, Z_adj_TS])
X <- as.matrix(LassoUKbiobank[, -Z_adj_TS])


set.seed(21697626)
train <- sample(1:nrow(LassoUKbiobank), 0.7 * nrow(LassoUKbiobank))
test <- seq(1, nrow(LassoUKbiobank))[-train]

set.seed(100)
model_lasso <- cv.glmnet(X[train, ], Y[train], alpha = 1,
                         family = "gaussian")
plot(model_lasso)

model_lasso$lambda.min

model_lasso$lambda.1se

min(model_lasso$cvm)
  