## Plotting the p-values of the missing telomeres:
rm(list=ls())


webshot::install_phantomjs()

library(MASS) # to access Animals data sets
library(scales)
library(ggplot2)
library(dplyr)
library(kableExtra)

# load the files containing the p-vals of the telomeres::
setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outputs/Results")
for (i in 1:7){
  if(i==1){
    
    Pvalues<-readRDS(paste0("univ_pvalues_several_nodes_one_core_", i, ".rds"))
  }
  else{
    Pvalues<-rbind(Pvalues,readRDS(paste0("univ_pvalues_several_nodes_one_core_", i, ".rds")))
    }
}

#length(data$Sex)-summary(model1)$df[2]-summary(model1)$df[1]

#c(rownames(low),rownames(mid))

plot_betas<-Pvalues %>% select(-beta) %>% mutate( p_val=ifelse (p_val==0,8.103441e-300,p_val))

plot_betas$number = 1:length(plot_betas[,1])

annotation <- data.frame(x = c(15,15),
                         y = c(20,0),
                         label = c("bonferroni correction", "0.05 threshold"))

pdf("rplot_p_vals.pdf")

ggplot(plot_betas,aes(y=-log(p_val),number,col=-log(p_val)>-log(0.05/length(plot_betas[,1])))) +
  scale_color_discrete(name = "p-vals", labels = c("not significant", "significant")) +
  geom_label(data =annotation, aes( x=x, y=y, label=label),color="black", size=3) +
  geom_point() + 
  geom_hline(yintercept=-log(0.05/length(plot_betas[,1])), linetype="dashed", color = "blue") +
  geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red") +
  scale_y_sqrt() +
  xlab("variables")+ylab("p values (negative Log) ")+labs(title="P-values of univariate regression")
dev.off()   

Pvalues$number = 1:length(Pvalues[,1])

pdf("rplot_beta.pdf")
P_vals<-Pvalues%>%filter(beta<4)
ggplot(Pvalues,aes(y=-log(p_val),x=beta,col=-log(p_val)>-log(0.05/length(plot_betas[,1])))) +
  scale_color_discrete(name = "p-vals", labels = c("not significant", "significant")) +
  geom_point()+xlab("coefficient")+ylab("p values (negative Log) ")+labs(title="Coefficient vs p-value")

ggplot(P_vals,aes(y=-log(p_val),x=beta,col=-log(p_val)>-log(0.05/length(plot_betas[,1])))) +
  scale_color_discrete(name = "p-vals", labels = c("not significant", "significant")) +
  geom_point()+xlab("coefficient")+ylab("p values (negative Log) ")+
  labs(title="Excluding sex, Coefficient vs p-value")
dev.off()   

#Identify those with significant p_vals:

Pvalues %>% filter(p_val<0.05/length(plot_betas[,1])) %>% 
            arrange(desc(abs(beta))) %>% 
            select(var_names,beta) %>% 
            filter(var_names !="Z_adj_TS.1") %>% kable(.,"html") %>%
  kable_styling("striped", full_width = FALSE, htmltable_class = 'lightable-classic-2') %>%
  kable_paper() %>%
save_kable("rplot_table.png", bs_theme = "flatly")
