library(tidyverse)
library(focus)
library(withr)
library(gridExtra)

ICU_ukbiobank<-readRDS("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/ANALYSIS_STATISTICS/univariate_covid_ukbb/Advanced_Analysis/icu_biobank.rds")
#Sensitivity Analysis - use Lasso analysis and stability
confint(lm(Leukocyte_telomere_length~leukocyte_count,ICU_ukbiobank))

#Sensitivity Analysis - use Lasso analysis and stability
colnames(ICU_ukbiobank)
#deaths
y<- ICU_ukbiobank%>% filter(death==1) %>% select(case) %>% mutate(case=as.numeric(case)-1)
x<- ICU_ukbiobank%>% filter(death==1) %>% select(-case,-case,-eid) 
# %>% select(Age,Leukocyte_telomere_length,Sex)

x_train <- model.matrix( ~ .-1, x)
glm<-cv.glmnet(x = x_train, y = as.matrix(y), penalty.factor=c(rep(1,ncol(x_train))), family = "binomial")

pdf(file = paste0())
plot(glm)
dev.off()

beta_lasso=coef(glm, s="lambda.min")[2:(ncol(x_train)+1),]
selected_lasso=names(beta_lasso)[which(beta_lasso!=0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

out = VariableSelection(x = x_train,seed=1,y = as.matrix(y),family = "binomial")

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

par(mar=c(10,5,1,1))

selection_plot<-data.frame(names=names(selprop),prop=selprop)%>%mutate(good=as.factor(ifelse(selprop>hat_params[2],"1","0"))) %>% arrange(desc(prop))

a <- ifelse(selection_plot$good == "0", "#999999", "#3c7be8")

ggplot(selection_plot,aes(x=fct_reorder(names,desc(prop)),y=prop,col=good,fill=good))+
  geom_segment(aes(yend=0,xend=names),size=2,lineend = "round")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1,colour=a),
        text = element_text(size=13),
        strip.placement = "outside",
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))


#############################################################################################################
# ICU

#Sensitivity Analysis - use Lasso analysis and stability
variables_icu<-c("bressupdays","aressupdays",
                 "bcardsupdays","acardsupdays","rensupdays",
                 "neurosupdays","gisupdays","dermsupdays",
                 "orgsupmax","bsc_res_sup","adv_res_sup",
                 "other_organ_support","other_organ_support")
#icu
y<- ICU_ukbiobank%>% filter(icu_stay==1) %>% select(case) %>% mutate(case=as.numeric(case)-1)
x<- ICU_ukbiobank%>% filter(icu_stay==1) %>% select(-variables_icu)%>% select(-case,-eid) 

# %>% select(Age,Leukocyte_telomere_length,Sex)

x_train <- model.matrix( ~ .-1, x)
glm<-cv.glmnet(x = x_train, y = as.matrix(y), penalty.factor=c(rep(1,ncol(x_train))), family = "binomial")

pdf(file = paste0())
plot(glm)
dev.off()

beta_lasso=coef(glm, s="lambda.min")[2:(ncol(x_train)+1),]
selected_lasso=names(beta_lasso)[which(beta_lasso!=0)]
print(paste0(length(selected_lasso), " proteins are selected"))
print(selected_lasso)

out = VariableSelection(x = x_train,seed=1,y = as.matrix(y),family = "binomial")

pdf()
CalibrationPlot(out)
dev.off()

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

par(mar=c(10,5,1,1))

selection_plot<-data.frame(names=names(selprop),prop=selprop)%>%mutate(good=as.factor(ifelse(selprop>hat_params[2],"1","0"))) %>% arrange(desc(prop))

a <- ifelse(selection_plot$good == "0", "#999999", "#3c7be8")

ggplot(selection_plot,aes(x=fct_reorder(names,desc(prop)),y=prop,col=good,fill=good))+
  geom_segment(aes(yend=0,xend=names),size=2,lineend = "round")+
  geom_hline(aes(yintercept = hat_params[2]),color="#990000",linetype="dashed")+
  scale_fill_manual(values=c("#999999", "#56B4E9"))+
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  ylab("Selection proportion")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1,colour=a),
        text = element_text(size=13),
        strip.placement = "outside",
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))



########################################################################################################################################################

# Age split

ICU_ukbiobank<-ICU_ukbiobank %>% mutate( Age_50_60 = as.factor(ifelse(Age>50 & Age <=60,"1","0")),
                                         Age_60_70 = as.factor(ifelse(Age>60 & Age <=70,"1","0")),
                                         Age_70_over = as.factor(ifelse(Age>70,"1","0")))

# ICU regressions:
summ_icu<-summary(glm(icu_stay~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
              filter(Age_50_60==1)%>% select(-case,-death)))
summ_icu_1<-data.frame(coef(summ_icu)) %>% mutate(names=rownames(coef(summ_icu)),group="Age_50_60")

#Death regressions:

summary(glm(death~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=binomial,data=ICU_ukbiobank%>%
              filter(Age_50_60==1) %>% select(-case)))
summary(glm(case~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
              filter(Age_50_60==1) %>% filter(death==1)))

summ_death<-summary(glm(case~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                          filter(Age_50_60==1) %>%filter(death==1)))
summ_death_1<- data.frame(coef(summ_death)) %>% mutate(names=rownames(coef(summ_death)),group="Age_50_60")

########################################################################################################################################################

# ICU regressions:
summ_icu<-summary(glm(icu_stay~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                        filter(Age_60_70==1)%>% select(-case,-death)))
summ_icu_2<-data.frame(coef(summ_icu)) %>% mutate(names=rownames(coef(summ_icu)),group="Age_60_70")

#Death regressions:

summ_death<-summary(glm(case~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                          filter(Age_60_70==1) %>%filter(death==1)))
summ_death_2<- data.frame(coef(summ_death)) %>% mutate(names=rownames(coef(summ_death)),group="Age_60_70")


pd_2<-ggplot(summ_death,aes(x=names,y=exp(Estimate)))+
  geom_point() +
  ggtitle("Death relative Ukbiobank and Covid cases")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), width=0)+
  xlab("Variables and levels") +
  ylab("OR for relative death vs death in the biobank")

pi_2<-ggplot(summ_icu,aes(x=names,y=exp(Estimate)))+
  xlab("Variables and levels") +
  ylab("OR for relative ICU admissions vs historical admissions in ICU in the biobank")+
  ggtitle("Icu Admissions")+
  geom_point() + 
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), width=0)

########################################################################################################################################################

# ICU regressions:
summ_icu<-summary(glm(icu_stay~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                        filter(Age_70_over==1)%>% select(-case,-death)))
summ_icu_3<-data.frame(coef(summ_icu)) %>% mutate(names=rownames(coef(summ_icu)),group="Age_70_over")

#Death regressions:

summ_death<-summary(glm(case~Sex+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                          filter(Age_70_over==1) %>%filter(death==1)))
summ_death_3<- data.frame(coef(summ_death)) %>% mutate(names=rownames(coef(summ_death)),group="Age_70_over")

summ_death <- rbind(summ_death_1,summ_death_2,summ_death_3)

pd <- position_dodge(0.5) 

ggplot(summ_death %>%filter(names!="(Intercept)"), aes(x=names, y=exp(Estimate), colour=group, group=group)) + 
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), colour="black", width=0.3, position=pd)+
  geom_point(position=pd, size=1, shape=21, fill="white") + # 21 is filled circle
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  xlab("Variables and levels") +
  ylab("OR for relative death vs death in the biobank")+
  scale_colour_hue(name="Age group",    # Legend label, use darker colors
                   breaks=c("Age_50_60", "Age_60_70","Age_70_over"),
                   labels=c("[50,60)", "[60,70)",">70"),
                   l=20) +                    # Use darker colors, lightness=40
  ggtitle("Death relative Ukbiobank and Covid cases")+
  expand_limits(y=0) +                        # Expand y range
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.95,0.76),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic")
  ) # Position legend in bottom right

### ICU plots

summ_icu<-rbind(summ_icu_1,summ_icu_2,summ_icu_3)

ggplot(summ_icu %>%filter(names!="(Intercept)"), aes(x=names, y=exp(Estimate), colour=group, group=group)) + 
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), colour="black", width=0.3, position=pd)+
  geom_point(position=pd, size=1, shape=21, fill="white") + # 21 is filled circle
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  xlab("Variables and levels") +
  ylab("OR for relative icu vs icu cases in the biobank")+
  scale_colour_hue(name="Age group",    # Legend label, use darker colors
                   breaks=c("Age_50_60", "Age_60_70","Age_70_over"),
                   labels=c("[50,60)", "[60,70)",">70"),
                   l=20) +                    # Use darker colors, lightness=40
  ggtitle("ICU admissions relative Ukbiobank and Covid cases")+
  expand_limits(y=0) +                        # Expand y range
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.95,0.76),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic")
  ) # Position legend in bottom right

########################################################################################################################################################

# Sex split - Select Women
summ_icu_F<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                        filter(Sex=="Female")%>% select(-death)))
summ_icu_M<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                          filter(Sex=="Male")%>% select(-death)))


summ_death_F<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                            filter(Sex=="Female") %>%filter(death==1)))
summ_death_M<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                            filter(Sex=="Male") %>%filter(death==1)))

summ_death_F<- data.frame(coef(summ_death_F)) %>% mutate(names=rownames(coef(summ_death_F)))
summ_death_M<- data.frame(coef(summ_death_M)) %>% mutate(names=rownames(coef(summ_death_M)))


#####################################

summ_death_F<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                            filter(Sex=="Female") %>%filter(death==1)))

summ_death<- data.frame(coef(summ_death)) %>% mutate(names=rownames(coef(summ_death)))

summ_death_M<-summary(glm(case~Age+Ethnicity+Diabetes+leukocyte_count+Average_household_income+Vascular_disease+Smoking_status+Alcohol_drinker_status+Leukocyte_telomere_length,family=quasibinomial,data=ICU_ukbiobank%>% 
                            filter(Sex=="Male") %>%filter(death==1)))

summ_death_2<- data.frame(coef(summ_death_M)) %>% mutate(names=rownames(coef(summ_death_M)),group="Male")
summ_death_1<- data.frame(coef(summ_death_F)) %>% mutate(names=rownames(coef(summ_death_F)),group="Female")

sum_death_sex<-rbind(summ_death_1,summ_death_2)

pd <- position_dodge(0.5) 

ggplot(sum_death_sex %>%filter(names!="(Intercept)"), aes(x=names, y=exp(Estimate), colour=group, group=group)) + 
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), colour="black", width=0.3, position=pd)+
  geom_point(position=pd, size=1, shape=21, fill="white") + # 21 is filled circle
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  xlab("Variables and levels") +
  ylab("OR for relative death vs death in the biobank")+
  scale_colour_hue(name="Sex",    # Legend label, use darker colors
                   breaks=c( "Male","Female"),
                   labels=c("Male", "Female"),
                   l=20) +                    # Use darker colors, lightness=40
  ggtitle("Death relative UK-biobank and Covid cases")+
  expand_limits(y=0) +                        # Expand y range
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.955,0.81),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic")
  ) # Position legend in bottom right

##For ICU

summ_icu_M<- data.frame(coef(summ_icu_M)) %>% mutate(names=rownames(coef(summ_icu_M)),group="Male")
summ_icu_F<- data.frame(coef(summ_icu_F)) %>% mutate(names=rownames(coef(summ_icu_F)),group="Female")

sum_icu_sex<-rbind(summ_icu_M,summ_icu_F)


ggplot(sum_icu_sex %>%filter(names!="(Intercept)"), aes(x=names, y=exp(Estimate), colour=group, group=group)) + 
  geom_errorbar(aes(ymin=exp(Estimate-1.95*abs(Std..Error)), ymax=exp(Estimate+1.95*abs(Std..Error))), colour="black", width=0.3, position=pd)+
  geom_point(position=pd, size=1, shape=21, fill="white") + # 21 is filled circle
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  xlab("Variables and levels") +
  ylab("OR for relative icu vs icu cases in the biobank")+
  scale_colour_hue(name="Sex",    # Legend label, use darker colors
                   breaks=c( "Male","Female"),
                   labels=c("Male", "Female"),
                   l=20) +                    # Use darker colors, lightness=40
  ggtitle("ICU admissions relative UK-biobank and Covid cases")+
  expand_limits(y=0) +                        # Expand y range
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.955,0.81),
        text = element_text(size=10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size=14, face="bold.italic")
  ) # Position legend in bottom right

