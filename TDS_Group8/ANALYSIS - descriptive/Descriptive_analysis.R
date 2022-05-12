setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording")
UKbiobank<- readRDS("AA_Group_data/DATASET/data_ukbiobank.rds")
install.packages("table1")
require(table1)
table1(~ Z_adj_TS + Sex + Age_recruited + BMI + Location + avg_household_income + CurrSmok + IMD + disability_allowance + Age_HBP + heart_attack_bin, data=UKbiobank)

library(tidyverse)
UKbiobank %>% ggplot(aes(x = Z_adj_TS)) + geom_histogram(bins=100, color="darkblue", fill="lightblue")

