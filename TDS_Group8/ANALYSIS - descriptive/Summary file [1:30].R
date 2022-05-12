args <- commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

library(ggplot2)
library(dplyr)
# library(naniar)


setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outputs/plots")
data<-readRDS("../official_dataset.rds")
annot<-readRDS("../annot.rds")

#Remove those which have Telomere length missing, unless there's a decent amount, then impute??



data<- data %>% mutate(
  IMD = ifelse(is.na(IMD_england),
               ifelse(is.na(IMD_scotland),IMD_wales,IMD_scotland),
               IMD_england),
  Location=ifelse(is.na(IMD_england),
                  ifelse(is.na(IMD_scotland),"Wales","Scotland"),
                  "England"),
  birth_weight = ifelse(birth_weight>5,birth_weight/2.205,birth_weight),
  distance_home_work = ifelse(distance_home_work>=1000,NA,distance_home_work),
  dur_mod_exc = ifelse(dur_mod_exc/1440>=0.75,NA,dur_mod_exc),
  dur_vig_exc = ifelse(dus_vig_exc/1440>=0.75,NA,dus_vig_exc),
  SleepDur = ifelse(SleepDur>15,NA,SleepDur),
  YOB = 2021 - YOB
)


data_noTS_NA <- data %>% filter(is.na(Z_adj_TS)!=TRUE) %>% select( - c(IMD_scotland,IMD_wales,IMD_scotland)) %>% mutate(
  birth_weight = ifelse(birth_weight>5,birth_weight/2.205,birth_weight))

summary(data$Z_adj_TS)
sum(is.na(data$Z_adj_TS))/length(data[,1])

data <- data %>% select( - c(IMD_scotland,IMD_wales,IMD_scotland))
col_num<-match("Z_adj_TS",colnames(data))
subset<-cbind(data[,1:30],data[col_num])


##########################
# PLOTTING
colnames(data[,1:39])
#Sex
plot(data[1])
#Age
hist(2021-data[,2])
#IMD_town
mean(data[,3],na.rm=TRUE)
hist(data[,3]) #the data is slightly negatively skewed

#IMD, and by areas
hist(data[,4])

ggplot(data,aes(x=IMD,col=Location))+geom_density(position = "identity", alpha = 0.2, bins = 50)
plot(data$avg_household_income)   #See it's not numeric & many not answers
hist(data$working_week)

#Age recruited
hist(data$Age_recruited)

#birthweight
max(data$birth_weight,na.rm=TRUE)
hist(data$birth_weight)

#most likely these are outliers and should be transformed into kg from pounds
sort(data$birth_weight[is.na(data$birth_weight)==FALSE],decreasing=TRUE)
sort(data$birth_weight[is.na(data$birth_weight)==FALSE],decreasing=FALSE)

hist(data$birth_weight)


hist(data$birth_weight)
# number of cancers
hist(data$num_cancers)

# Acc lived & own accomodation / length @address - use with pollution & IMD
plot(data$accomodation_lived)

plot(data$own_accomodation) #should relabel this?
hist(data$length_current_address)

# distance from work/nightshits/length + exercice
levels(data$Night_shift)
plot(data$Night_shift)
hist(data$working_week)

sort(data$working_week[is.na(data$working_week)==FALSE],decreasing=TRUE)/7

max(data$distance_home_work,na.rm = TRUE)
sort(data$distance_home_work[is.na(data$distance_home_work)==FALSE],decreasing=TRUE)
hist(data$distance_home_work) #Need to remove outliers?

hist(data$num_days_walk_10mins)

hist(data$num_days_mod_exc)
hist(data$dur_mod_exc)

hist(data$num_says_vig_exc)
hist(data$dus_vig_exc)

plot(data$walk_pace)

#Sleeping:
sort(data$SleepDur,decreasing=TRUE)

hist(data$SleepDur)
levels(data$narcolepsy)


plot(data$narcolepsy) + title("narcolepsy")
plot(data$Insomnia) + title("Insomnia")

#Smoking:
plot(data$CurrSmok)
plot(data$PastSmok)
plot(data$n_smoke_hshld)
hist(data$tobacco_exp_hshld)
hist(data$tobacco_exp_outside_hshld)

# Food_stuff:
colnames(data)[31:35]

hist(data$cooked_veggie_intake)
hist(data$raw_veggie_intake)
hist(data$fresh_fruit_intake)
hist(data$dried_fruit_intake)
plot(data$oily_fish_intake)

###########################
#Available data: 
length(data$Sex)-apply(data,MARGIN=2,function(x){sum(is.na(x))})


###########################
# Analysis of the betas of influence:

max_levels<-1
for (i in 1:length(data)){
  max_levels<-  ifelse(length(levels(data[,i]))<max_levels,max_levels,length(levels(data[,i])))
}

Beta = Pvalues = Names= matrix(NA, nrow = max_levels,ncol = (length(data)-1))

## split this into different nodes and run better

data<-data %>% select(-c("lymp_count","lymp_count_dev","mon_count","mon_count_dev","lymp_dev",
                         "myocardial_date","cancer_date","stdate_dementia","Unadj_TS","Adj_TS"))

ids=as.character(cut(1:length(data), breaks = nchunks, labels = 1:nchunks))

pvalues=data.frame(cbind(data[,ids == ichunk ],Z_adj_TS=data$Z_adj_TS))

colnames(pvalues)[1:length(pvalues)]


for (j in  1:length(pvalues)) {
print(j)
    if(class(pvalues[,j])=="factor" & colnames(pvalues[j])!="Z_adj_TS"){
  
  len<-length(levels(pvalues[,j]))

  model1 = lm(Z_adj_TS~pvalues[,j],data=pvalues)
  Beta[, j] = c(coefficients(model1)[2:c(len)],rep(NA,max_levels-len+1))
  Names[, j] = c(row.names(summary(model1)$coefficients)[2:c(len)],rep(NA,max_levels-len+1))
  Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                   rep(NA,max_levels-c(length(summary(model1)$coefficients[,1]))+1))
    }
    if(colnames(pvalues[j])!="Z_adj_TS"){
    model1 = lm(Z_adj_TS~pvalues[,j],data=pvalues)
    Beta[, j] = c(coefficients(model1)[2],rep(NA,max_levels-1))
    Names[, j] = c(colnames(pvalues[j]),rep(NA,max_levels-1))
    Pvalues[, j] = c(summary(model1)$coefficients[2:c(length(summary(model1)$coefficients[,1])),4],
                     rep(NA,max_levels-c(length(summary(model1)$coefficients[,1]))+1))
}}

plot_betas<-data.frame("var_names"=c(Names),"p_val"=as.numeric(Pvalues),"beta"=c(Beta)) %>% filter(is.na(var_names)==FALSE)

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS")
ifelse(dir.exists("Results"),"",dir.create("Results"))
saveRDS(plot_betas, paste0("Results/univ_pvalues_several_nodes_one_core_", ichunk, ".rds"))


########################

# See % of missingness in these variables

# % distribution and variable type

# See if variables can be aggregated or unknown patters?

# Qualify the variables and see which ones might be the most interesting - NA analysis from these

#Decide imputation/removal/ignoring - MNAR/MCAR/other

Data<-data[1:1000,]
na_perc<-round(100*sapply(data, function(x) sum(length(which(is.na(x)))))/length(data$Sex),3)
rownames(data.frame(na_perc))

#Good variables
data.frame(na_perc[na_perc<25])

#Middling variables
na_perc[na_perc<50 &na_perc>25]

#poor Quality variables
mid<-data.frame(na_perc[na_perc>50 & na_perc<=70])
low<-data.frame(na_perc[na_perc>70])

filter_missings <- Data %>% select(-c(rownames(low),rownames(mid)))

levels(data$menopause_bi)

########################


