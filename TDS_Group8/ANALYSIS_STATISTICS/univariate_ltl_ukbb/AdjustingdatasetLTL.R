setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")


dff3<- readRDS("DATASET/data_ukbiobank.rds")

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_ltl_ukbb")

saveRDS(dff3, "tempdataset.rds")

dff3<- readRDS("tempdataset.rds")


#  edit variable names

dff3 <- rename(dff3, Leukocyte_telomere_length = Z_adj_TS,
               Index_of_multiple_deprivation = IMD,
               Number_of_cancers = num_cancers,
               Accomodation_type = accomodation_lived,
               Narcolepsy = narcolepsy,
               Oily_fish_intake = oily_fish_intake,
               Non_oily_fish_intake = non_oily_fish_intake,
               Processed_meat_consumption = ProcMeat,
               Salt_consumption = salt_in_food,
               Coffee_type = coffee_type_intake,
               Long_standing_illness = long_standing_illness,
               Disability_allowance = disability_allowance,
               Myocardial_disease_diagnosis_type = souce_myocard,
               Accomodation_ownership = own_accomodation,
               Tea_intake = tea_intake,
               Raw_veggie_intake = raw_veggie_intake,
               Num_days_walk_10mins = num_days_walk_10mins,
               Had_cancer = cancer_bin,
               Cancer_diagnosis_bydoctor = cancer_doctor,
               Fresh_fruit_intake = fresh_fruit_intake,
               Length_at_current_address = length_current_address,
               Average_household_income = avg_household_income,
               Days_per_wk_moderate_activity = num_days_mod_exc,
               Days_per_wk_vigorous_activity = num_says_vig_exc,
               Walk_pace = walk_pace,
               Smoking_status = CurrSmok,
               Smokers_in_household = n_smoke_hshld,
               Bread_consumption = breadcons,
               Mother_smoking_birth = mat_smoking_birth,
               Diabetes = diabetes_doctor,
               Menopause = menopause_bi,
               Contraceptive_pill = contra_pill,
               Hormone_replacement_therapy = hormone_rep_thera,
               Age_highBP_diagnosed = Age_HBP,
               Stillbirths = num_stillbirths,
               Miscarriages = num_miscarriages,
               Num_vitamin_supplements = vit_suplements,
               Illnesses_of_father = illnesses_father,
               Illnesses_of_mother = illnesses_mother,
               Illnesses_of_siblings = illnesses_siblings,
               Alcohol_drinker_status = alcohol_drink_status,
               Age_last_took_cannabis = LastWeedAge,
               No_eggsdairywheatorsugar = NoEDWS,
               Diagnosed_with_common_disease = common_disease_diag,
               Level_of_NO2_air_pollution = NO2_poll, 
               Diastolic_blood_pressure = Dia_BP_auto,
               Systolic_blood_pressure = Syt_BP_auto,
               Number_in_household = numb_household,
               Heart_attack = heart_attack_bin,
               Medic_comorbidities = medic_comorbidities_bin,
               High_blood_pressure = Age_HBP_bin,
               Myocardial_disease = myocardial_binary,
               Dementia = dementia_binary,
               Depressed_in_last2wks = depressed_last_2wks)

# relevel so reference is the "norm"/ without condition 


dff4 <- dff3 %>%mutate(
  Accomodation_ownership=relevel(Accomodation_ownership,ref="Rent"),
    Narcolepsy=relevel(Narcolepsy,ref="Never/rarely"),
    Diabetes=relevel(Diabetes,ref="No"),
  Cancer_diagnosis_bydoctor=relevel(Cancer_diagnosis_bydoctor,ref="No"),
  Menopause=relevel(Menopause,ref="No"),
  Contraceptive_pill=relevel(Contraceptive_pill,ref="No"),
  Hormone_replacement_therapy=relevel(Hormone_replacement_therapy,ref="No"),
    Stillbirths=relevel(Stillbirths,ref="None"),
  Myocardial_disease_diagnosis_type=relevel(Myocardial_disease_diagnosis_type,ref="Not diagnosed"),
  Depressed_in_last2wks=relevel(Depressed_in_last2wks,ref="Not at all"),
Walk_pace=relevel(Walk_pace, ref = "Steady average pace"),
Salt_consumption=relevel(Salt_consumption, ref = "Sometimes"),
Mother_smoking_birth=relevel(Mother_smoking_birth, ref = "No"),
Illnesses_of_father = relevel(Illnesses_of_father, ref = "None of the above"),
Illnesses_of_mother = relevel(Illnesses_of_mother, ref = "None of the above"), 
Illnesses_of_siblings = relevel(Illnesses_of_siblings, ref = "None of the above")) 

#   Change factors with 0 / 1 to yes / no

dff5 <-dff4 %>%mutate(
  Had_cancer = as.factor(ifelse(Had_cancer == 0,"No","Yes")),

  Medic_comorbidities = as.factor(ifelse(Medic_comorbidities == 0,"No","Yes")),
  
  High_blood_pressure = as.factor(ifelse(High_blood_pressure == 0,"No","Yes")),
  
  Heart_attack = as.factor(ifelse(Heart_attack == 0,"No","Yes")),
  
  Myocardial_disease = as.factor(ifelse(Myocardial_disease == 0,"No","Yes")),
  
  Dementia = as.factor(ifelse(Dementia == 0,"No","Yes")))



#   dia_BP (also an auto), age diabetes bin and menopause age bin (repeats of same info)

dff6 <- subset(dff5, select = -c(Age_diabetes_bin, MenopauseAge_bin, source_dementia, Syst_BP, Dia_BP))


saveRDS(data.frame(dff6), "tempdataset.rds")



#Updated Variable groups to new names: 

Lifestyle<-c("Num_days_walk_10mins","Days_per_wk_moderate_activity","Days_per_wk_vigorous_activity","Walk_pace","SleepDur",
             "Smoking_status","PastSmok","Smokers_in_household","tobacco_exp_hshld","Location","Sunburn_childhood","Mother_smoking_birth","Alcohol_drinker_status","Level_of_NO2_air_pollution","PM10_poll", "Age_last_took_cannabis")
Lifestyle<-data.frame(variables=Lifestyle,type=rep("Lifestyle",length(Lifestyle)))

Dietary<-c("cooked_veggie_intake","Raw_veggie_intake","Fresh_fruit_intake",
           "dried_fruit_intake","Oily_fish_intake","Non_oily_fish_intake","Processed_meat_consumption","cheese_cons","Bread_consumption",
           "Salt_consumption","Tea_intake","Coffee_type","Water_intake","Alcohol_intake_freq","Num_vitamin_supplements", "No_eggsdairywheatorsugar")
Dietary<-data.frame(variables=Dietary,type=rep("Dietary",length(Dietary)))

Socio_economic<-c("IMD","Index_of_multiple_deprivation", "Accomodation_type","Accomodation_ownership","Length_at_current_address","Average_household_income",
                  "Disability_allowance","Number_in_household")
Socio_economic<-data.frame(variables=Socio_economic,type=rep("Socio_economic",length(Socio_economic)))

Medical<-c("leukocyte_count","Number_of_cancers","Long_standing_illness","Insomnia","Narcolepsy", "Had_cancer", "Diagnosed_with_common_disease",
           "Diabetes","Cancer_diagnosis_bydoctor","presc_medicine","Menopause","Contraceptive_pill","Hormone_replacement_therapy","Age_highBP_diagnosed",
           "forces_vital_capacity","MenopauseAge","Stillbirths","Miscarriages","Illnesses_of_father","Illnesses_of_mother",
           "Illnesses_of_siblings","Never_eat_eggs_dairy_wheat_sugar","common_disease_diag","Myocardial_disease_diagnosis_type",
           "Dementia","Systolic_blood_pressure","Diastolic_blood_pressure","Heart_attack","Medic_comorbidities","Age_diabetes_bin",
           "MenopauseAge_bin","High_blood_pressure","Myocardial_disease","Depressed_in_last2wks","Cystatin_C","Vascular_disease", "Leukocyte_telomere_length")

Medical<-data.frame(variables=Medical,type=rep("Medical",length(Medical)))


Personal<-c("Sex","Age","Ethnicity","BMI","Height")
Personal<-data.frame(variables=Personal,type=rep("Personal",length(Personal)))


subgroup_organization_updatedLTL<-left_join(data.frame(variables=colnames(dff6)),rbind(Lifestyle,Dietary,Socio_economic,Medical,Personal))
saveRDS(subgroup_organization_updatedLTL, "subgroup_organization_updatedLTL.rds")

