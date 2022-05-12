#
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

##Added a line in the regression function to print the standard deviation
# also modified the chunks since can be run without the submit
# saved result as well without the number of the chunk
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)


## DATA LOADING



setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")



dff3<- readRDS("DATASET/data_ukbiobank_updated.rds")  

dff3a<-subset(dff3, Sex!="Male")  %>% select(-Sex)



###########
### LTL

#break it down so it's easier to run:


#See "covid_plot_univar.sh" to see how to submit the file with the chunks data

#break into chunks or silence if necessary:

dff4<-dff3a %>% select(-Leukocyte_telomere_length) %>% select(-eid) %>% select(-Age)
ids=as.character(cut(1:length(dff4)-2, breaks = nchunks, labels = 1:nchunks))

#dff4[,ids == ichunk] <- sub this back in the dataframe below when submiting is back on
pvalues=data.frame(dff4 , "Leukocyte_telomere_length" = dff3a %>% 
                     select(Leukocyte_telomere_length) %>% 
                     mutate(Leukocyte_telomere_length = as.numeric(Leukocyte_telomere_length)),
                   "Age" = dff3a %>% 
                     select(Age) %>% 
                     mutate(Age = as.numeric(Age)))
                  
#divide factors and continuous:
Xvar_cts<-data.frame(pvalues[,!sapply(pvalues,is.factor)] %>% select(-Leukocyte_telomere_length, -Age))

Xvar_factors<- data.frame(pvalues[,sapply(pvalues,is.factor)])


#Check each variable one by one to see which significant:

Betas_pvals<-function(X){
  model1 = summary(lm(Leukocyte_telomere_length ~ X + Age, data=pvalues))$coefficients
  Beta = model1[2:nrow(model1),"Estimate"]
  Pvalues = model1[2:nrow(model1),"Pr(>|t|)"]
  varbl_names = rownames(model1)[2:nrow(model1)]
  sd = model1[2:nrow(model1),"Std. Error"]
  df<-data.frame("Variables"=varbl_names,"beta"=Beta,"pvals"=Pvalues,"sd"=sd)
  return(df)
}


## Generalized ratio test -Extension, if happy remove the # and run
lrt__a<- function(Xvar) {
  mod<-summary(lm(Z_adj_TS ~ Xvar, data=pvalues))
  LRT_ltl<- pchisq(mod$null-mod$deviance, mod$df.null-mod$df.residual,lower.tail = F)
  return(LRT_ltl)
}

## Nels work

if (length(Xvar_cts)!=0) {cts<-apply(Xvar_cts,2,FUN=Betas_pvals)}

if (length(Xvar_factors)!=0) {factors<-apply(Xvar_factors,2,FUN=Betas_pvals)}

ltl_results<-rbind(if (length(Xvar_factors)!=0) {do.call(rbind, Map(data.frame, factors))},
                   if (length(Xvar_cts)!=0) {do.call(rbind, Map(data.frame, cts))}) 


ltl_results<- ltl_results %>% mutate(Variables= rownames(ltl_results))

separated_ltl<-separate(ltl_results,Variables,into =c("Variable","factor"), sep=".X")


setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/AA_Group_data")

subgroup<-readRDS("./DATASET/ukbiobank_subgrouped.rds")

subgroup<-subgroup %>% filter(variables %in% colnames(pvalues)) %>% mutate(Variable= variables) %>% select(-variables) 

## Remove.Age.SexMale
separated_ltl<-left_join(separated_ltl,subgroup, by = "Variable") %>% mutate(factor=ifelse(is.na(factor),Variable,factor)) %>% 
  filter(!grepl('.Age', Variable)) 


saveRDS(data.frame(separated_ltl), paste0("../ANALYSIS_STATISTICS/univariate_ltl_ukbb/Results_univ_reg/univariate_ltl_adjusted_update_females.rds"))

debug_contr_error <- function (dat, subset_vec = NULL) {
  if (!is.null(subset_vec)) {
    ## step 0
    if (mode(subset_vec) == "logical") {
      if (length(subset_vec) != nrow(dat)) {
        stop("'logical' `subset_vec` provided but length does not match `nrow(dat)`")
      }
      subset_log_vec <- subset_vec
    } else if (mode(subset_vec) == "numeric") {
      ## check range
      ran <- range(subset_vec)
      if (ran[1] < 1 || ran[2] > nrow(dat)) {
        stop("'numeric' `subset_vec` provided but values are out of bound")
      } else {
        subset_log_vec <- logical(nrow(dat))
        subset_log_vec[as.integer(subset_vec)] <- TRUE
      } 
    } else {
      stop("`subset_vec` must be either 'logical' or 'numeric'")
    }
    dat <- base::subset(dat, subset = subset_log_vec)
  } else {
    ## step 1
    dat <- stats::na.omit(dat)
  }
  if (nrow(dat) == 0L) warning("no complete cases")
  ## step 2
  var_mode <- sapply(dat, mode)
  if (any(var_mode %in% c("complex", "raw"))) stop("complex or raw not allowed!")
  var_class <- sapply(dat, class)
  if (any(var_mode[var_class == "AsIs"] %in% c("logical", "character"))) {
    stop("matrix variables with 'AsIs' class must be 'numeric'")
  }
  ind1 <- which(var_mode %in% c("logical", "character"))
  dat[ind1] <- lapply(dat[ind1], as.factor)
  ## step 3
  fctr <- which(sapply(dat, is.factor))
  if (length(fctr) == 0L) warning("no factor variables to summary")
  ind2 <- if (length(ind1) > 0L) fctr[-ind1] else fctr
  dat[ind2] <- lapply(dat[ind2], base::droplevels.factor)
  ## step 4
  lev <- lapply(dat[fctr], base::levels.default)
  nl <- lengths(lev)
  ## return
  list(nlevels = nl, levels = lev)
}

debug_contr_error(pvalues)

pvalues$Bp_cholesterol_or_hormonal_medication