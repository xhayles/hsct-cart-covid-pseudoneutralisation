#load packages
require(ggplot2)
require(tidyverse)
require(lubridate)
library(gridExtra)
library(reshape2)
library(dplyr)
library(scales)
library(epitools)
library(writexl)
library(ggbreak)
library(graphics)
library(graphics)
library(survival)
library(plyr)
library(readr)
library(ggh4x)
library(ggpubr)
library(stats)
library(rstatix)
library(cowplot)

#load data
metadata <- read.csv("./data/dummy_metadata.csv")
nab <- read.csv("./data/dummy_nab_results.csv")
bindingab <- read.csv("./data/dummy_bindingab_results.csv")
elispot <- read.csv("./data/dummy_elispot_results.csv")

####nAb regression####
regression <- full_join(metadata,nab,by="pid")

regression$V1 <- as.Date(regression$V1)
regression$D0<- as.Date(regression$D0)
regression$visit_date <- as.Date(regression$ivig_date)
regression$ivig_date <- as.Date(regression$ivig_date)

regression$time_since_HSCT <- difftime(regression$V1, regression$D0, unit="days")

regression$time_since_HSCT2 <- ifelse(regression$time_since_HSCT <0, NA,
                                      ifelse(regression$time_since_HSCT <366, "<12 months",
                                             ifelse(regression$time_since_HSCT >365, ">12 months", "banana")))

regression$time_since_HSCT3 <- ifelse(regression$time_since_HSCT <183, "<6 months",
                                      ifelse(regression$time_since_HSCT <366, "<12 months",
                                             ifelse(regression$time_since_HSCT >365, ">12 months", "banana"))) 

regression$time_since_IVIG <- difftime(regression$visit_date, regression$ivig_date, unit="days")

regression_WT <- regression %>%
  filter(pseudovirus=="WT") #refers to "wild-type" i.e. B.1

regression_WT_postdose2 <- regression_WT %>%
  filter(visit_code=="postdose2")

regression_WT_postdose2_HC <- regression_WT_postdose2 %>%
  filter(serum_source == "PITCH") #create healthy control dataset for later

regression_WT_postdose2 <- regression_WT_postdose2 %>%
  filter(serum_source == "HSCT") #filter for HSCT recipients

regression_WT_postdose2$logic50 <- log10(regression_WT_postdose2$ic50)

regression_WT_postdose2$age2 <- ifelse(regression_WT_postdose2$age < 45, "<45",
                                       ifelse(regression_WT_postdose2$age < 65, "45-64", ">65"))

regression_WT_postdose2$neutralising <- as.factor(regression_WT_postdose2$neutralising)

regression_WT_postdose2$days_since_last_vac <- as.numeric(regression_WT_postdose2$days_since_last_vac)

regression_WT_postdose2$HSCT_group <- as.factor(regression_WT_postdose2$HSCT_group)
regression_WT_postdose2$HSCT_group <-relevel(regression_WT_postdose2$HSCT_group, ref="Auto")

regression_WT_postdose2$prime <- as.factor(regression_WT_postdose2$prime)
regression_WT_postdose2$prime <-relevel(regression_WT_postdose2$prime, ref="mRNA")

regression_WT_postdose2$serostatus <- as.factor(regression_WT_postdose2$serostatus)
regression_WT_postdose2$serostatus <-relevel(regression_WT_postdose2$serostatus, ref="N")

regression_WT_postdose2$indication <- as.factor(regression_WT_postdose2$indication)
regression_WT_postdose2$indication <-relevel(regression_WT_postdose2$indication, ref="Myeloid")

regression_WT_postdose2$prev_auto <- as.factor(regression_WT_postdose2$prev_auto)
regression_WT_postdose2$prev_auto <-relevel(regression_WT_postdose2$prev_auto, ref="Unchecked")

regression_WT_postdose2$prev_allo <- as.factor(regression_WT_postdose2$prev_allo)
regression_WT_postdose2$prev_allo <-relevel(regression_WT_postdose2$prev_allo, ref="Unchecked")

regression_WT_postdose2$no_prev_HSCT <- as.factor(regression_WT_postdose2$no_prev_HSCT)
regression_WT_postdose2$no_prev_HSCT <-relevel(regression_WT_postdose2$no_prev_HSCT, ref="Checked")

regression_WT_postdose2$lymph_count_vac2 <- as.numeric(regression_WT_postdose2$lymph_count_vac2)
regression_WT_postdose2$igg_level_vac2 <- as.numeric(regression_WT_postdose2$igg_level_vac2)

regression_WT_postdose2$lymph_depletion2 <- ifelse(regression_WT_postdose2$lymph_depletion == "None", "No",
                                                   ifelse(regression_WT_postdose2$lymph_depletion == "Anti-thymocyte globulin", "Yes",
                                                          ifelse(regression_WT_postdose2$lymph_depletion == "Campath", "Yes",
                                                                 ifelse(regression_WT_postdose2$lymph_depletion == "Post-transplant Cyclophosphamide", "Yes", 
                                                                        ifelse(regression_WT_postdose2$lymph_depletion == "", "Unknown",
                                                                               ifelse(regression_WT_postdose2$lymph_depletion == "Unknown", "Unknown",
                                                                                      regression_WT_postdose2$lymph_depletion))))))

regression_WT_postdose2$lymph_depletion3 <- ifelse(regression_WT_postdose2$lymph_depletion == "None", "No",
                                                   ifelse(regression_WT_postdose2$lymph_depletion == "Anti-thymocyte globulin", "ATG",
                                                          ifelse(regression_WT_postdose2$lymph_depletion == "Campath", "Non-ATG",
                                                                 ifelse(regression_WT_postdose2$lymph_depletion == "Post-transplant Cyclophosphamide", "Non-ATG", 
                                                                        ifelse(regression_WT_postdose2$lymph_depletion == "", "Unknown",
                                                                               ifelse(regression_WT_postdose2$lymph_depletion == "Unknown", "Unknown",
                                                                                      regression_WT_postdose2$lymph_depletion))))))


regression_WT_postdose2$lymph_depletion <- as.factor(regression_WT_postdose2$lymph_depletion)
regression_WT_postdose2$lymph_depletion <-relevel(regression_WT_postdose2$lymph_depletion, ref="None")

regression_WT_postdose2$lymph_depletion2 <- as.factor(regression_WT_postdose2$lymph_depletion2)
regression_WT_postdose2$lymph_depletion2 <-relevel(regression_WT_postdose2$lymph_depletion2, ref="No")

regression_WT_postdose2$lymph_depletion3 <- as.factor(regression_WT_postdose2$lymph_depletion3)
regression_WT_postdose2$lymph_depletion3 <-relevel(regression_WT_postdose2$lymph_depletion3, ref="Non-ATG")

regression_WT_postdose2$tbi <- as.factor(regression_WT_postdose2$tbi)
regression_WT_postdose2$tbi <-relevel(regression_WT_postdose2$tbi, ref="No")

regression_WT_postdose2$conditioning <-as.factor(regression_WT_postdose2$conditioning)
regression_WT_postdose2$conditioning <- relevel(regression_WT_postdose2$conditioning, ref="Reduced Intensity")

regression_WT_postdose2$ethnicity <- as.factor(regression_WT_postdose2$ethnicity)
regression_WT_postdose2$ethnicity <-relevel(regression_WT_postdose2$ethnicity, ref="White")

regression_WT_postdose2$diseasestatusattrialentry <- as.factor(regression_WT_postdose2$diseasestatusattrialentry)
regression_WT_postdose2$diseasestatusattrialentry <-relevel(regression_WT_postdose2$diseasestatusattrialentry, ref="Complete Remission")

regression_WT_postdose2$hsct_hadnolymphstmcllpriorvac <- as.factor(regression_WT_postdose2$hsct_hadnolymphstmcllpriorvac)
regression_WT_postdose2$hsct_hadnolymphstmcllpriorvac <-relevel(regression_WT_postdose2$hsct_hadnolymphstmcllpriorvac, ref="Checked")

regression_WT_postdose2$lymph_infusion_vac2 <- as.factor(regression_WT_postdose2$lymph_infusion_vac2)
regression_WT_postdose2$lymph_infusion_vac2 <- relevel(regression_WT_postdose2$lymph_infusion_vac2, ref="None")

regression_WT_postdose2$COVID_predose2 <- as.factor(regression_WT_postdose2$COVID_predose2)
regression_WT_postdose2$COVID_predose2 <- relevel(regression_WT_postdose2$COVID_predose2, ref="No")

regression_WT_postdose2$any_gvhd_vac1 <- as.factor(regression_WT_postdose2$any_gvhd_vac1)
regression_WT_postdose2$any_gvhd_vac1 <- relevel(regression_WT_postdose2$any_gvhd_vac1, ref="No/Unknown")

regression_WT_postdose2$any_gvhd_vac2 <- as.factor(regression_WT_postdose2$any_gvhd_vac2)
regression_WT_postdose2$any_gvhd_vac2 <- relevel(regression_WT_postdose2$any_gvhd_vac2, ref="No/Unknown")

regression_WT_postdose2$pred_5plus_vac1 <- as.factor(regression_WT_postdose2$pred_5plus_vac1)
regression_WT_postdose2$pred_5plus_vac1 <- relevel(regression_WT_postdose2$pred_5plus_vac1, ref="Unchecked") 

regression_WT_postdose2$pred_5plus_vac2 <- as.factor(regression_WT_postdose2$pred_5plus_vac2)
regression_WT_postdose2$pred_5plus_vac2 <- relevel(regression_WT_postdose2$pred_5plus_vac2, ref="Unchecked")

regression_WT_postdose2$time_since_HSCT2 <- as.factor(regression_WT_postdose2$time_since_HSCT2)
regression_WT_postdose2$time_since_HSCT2 <-relevel(regression_WT_postdose2$time_since_HSCT2, ref=">12 months")

regression_WT_postdose2$dosing_interval2 <- ifelse(regression_WT_postdose2$dosing_interval <42, "<6 weeks",
                                                   ifelse(regression_WT_postdose2$dosing_interval >41, ">6 weeks","banana"))   

####Univariate logistic regression model####

lapply(c("prime", "age", "age2","sex", "ethnicity", 
         "COVID_predose2",
         "time_since_HSCT2", "days_since_last_vac", 
         "dosing_interval2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","tbi",
         "lymph_depletion3", 
         "no_prev_HSCT",
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"
),
function(var) {
  
  
  formula = as.formula(paste("neutralising~", var))
  
  res.logit = glm(formula, family=binomial(link="logit"),  data=regression_WT_postdose2)
  
  summary(res.logit)
  
})

#ORs
lapply(c("prime", "age", "age2","sex", "ethnicity", 
         "COVID_predose2",
         "time_since_HSCT2", "days_since_last_vac", 
         "dosing_interval2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","tbi",
         "lymph_depletion3", 
         "no_prev_HSCT",
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"
),
function(var) {
  
  formula = as.formula(paste("neutralising~", var))
  
  res.logit = glm(formula, family=binomial(link="logit"), data=regression_WT_postdose2)
  
  exp(coef(res.logit))
  
})


#Get CIs for OR

lapply(c("prime", "age", "age2","sex", "ethnicity", 
         "COVID_predose2",
         "time_since_HSCT2", "days_since_last_vac", 
         "dosing_interval2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","tbi",
         "lymph_depletion3", 
         "no_prev_HSCT",
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"
),
function(var) {
  
  formula = as.formula(paste("neutralising~", var))
  
  res.logit = glm(formula, family=binomial(link="logit"), data=regression_WT_postdose2)
  
  exp(confint(res.logit))
  
})

###Multivariate LOGISTIC regression #####

multivariable_model_postdose2 = glm(neutralising~
                                      #factors with p <0.1
                                      lymph_count_vac2+
                                      age2+
                                      prime+
                                      HSCT_group+
                                      time_since_HSCT2+
                                      no_prev_HSCT+
                                      hastheparticipantreceivedrituxim+
                                      #a priori hypothesis
                                      dosing_interval2+
                                      COVID_predose2+
                                      indication,
                                      #conditioning+ #omit as data only available for allos
                                      #lymph_depletion3 #omit as data only available for allos
                                    family=binomial(link="logit"), data = regression_WT_postdose2)

summary(multivariable_model_postdose2)
exp(coef(multivariable_model_postdose2))
#exp(confint(multivariable_model_postdose2)) - CI too wide with dummy dataset

####alloHSCT ptna logistic regression####

regression_WT_postdose2_allos <- regression_WT_postdose2 %>%
  filter(HSCT_group == "Allo")
multivariable_model_postdose2_allos = glm(neutralising~
                                            #p <0.1
                                            age2+
                                            prime+
                                            time_since_HSCT2+
                                            no_prev_HSCT+
                                            hastheparticipantreceivedrituxim+
                                            lymph_count_vac2+
                                            #a priori hypothesis
                                            dosing_interval2+
                                            COVID_predose2+
                                            indication+
                                            conditioning+ #now included as data just for allos
                                            lymph_depletion3, #now included as data just for allos
                                          family=binomial(link="logit"), data = regression_WT_postdose2_allos)

summary(multivariable_model_postdose2_allos)
exp(coef(multivariable_model_postdose2_allos))
#exp(confint(multivariable_model_postdose2_allos)) - CI too wide with dummy dataset


###LOG10 LINEAR regression #####

summary(lm(logic50~prime, data=regression_WT_postdose2))
summary(lm(logic50~age, data=regression_WT_postdose2))
summary(lm(logic50~age2, data=regression_WT_postdose2))
summary(lm(logic50~sex, data=regression_WT_postdose2))
summary(lm(logic50~ethnicity, data=regression_WT_postdose2))
summary(lm(logic50~serostatus, data=regression_WT_postdose2))
summary(lm(logic50~time_since_HSCT, data=regression_WT_postdose2))
summary(lm(logic50~time_since_HSCT2, data=regression_WT_postdose2))
summary(lm(logic50~days_since_last_vac, data=regression_WT_postdose2))
summary(lm(logic50~HSCT_group, data=regression_WT_postdose2)) 
summary(lm(logic50~indication, data=regression_WT_postdose2))
summary(lm(logic50~diseasestatusattrialentry, data=regression_WT_postdose2))
summary(lm(logic50~conditioning, data=regression_WT_postdose2))
summary(lm(logic50~lymph_depletion, data=regression_WT_postdose2))
summary(lm(logic50~lymph_depletion2, data=regression_WT_postdose2))
summary(lm(logic50~lymph_depletion3, data=regression_WT_postdose2))
summary(lm(logic50~tbi, data=regression_WT_postdose2))
summary(lm(logic50~prev_auto, data=regression_WT_postdose2))
summary(lm(logic50~prev_allo, data=regression_WT_postdose2))
summary(lm(logic50~no_prev_HSCT, data=regression_WT_postdose2))
#summary(lm(logic50~toci_rec, data=regression_WT_postdose2)) 
summary(lm(logic50~hastheparticipantreceivedrituxim, data=regression_WT_postdose2))
summary(lm(logic50~lymph_count_vac2, data=regression_WT_postdose2))
summary(lm(logic50~igg_level_vac2, data=regression_WT_postdose2))
summary(lm(logic50~lymph_infusion_vac2, data=regression_WT_postdose2))
summary(lm(logic50~hsct_hadnolymphstmcllpriorvac, data=regression_WT_postdose2))
summary(lm(logic50~COVID_predose2, data=regression_WT_postdose2))
summary(lm(logic50~any_gvhd_vac1, data=regression_WT_postdose2))
summary(lm(logic50~any_gvhd_vac2, data=regression_WT_postdose2))
summary(lm(logic50~dosing_interval, data=regression_WT_postdose2))
summary(lm(logic50~dosing_interval2, data=regression_WT_postdose2))

multivariable_model_postdose2_loglinear = lm(formula = logic50~
                                               #p values <0.1
                                               lymph_count_vac2+
                                               age2+
                                               dosing_interval2+
                                               COVID_predose2+
                                               HSCT_group+
                                               time_since_HSCT2+
                                               indication+
                                               no_prev_HSCT+
                                               hastheparticipantreceivedrituxim+
                                               #a priori
                                               prime,
                                             data = regression_WT_postdose2)

summary(multivariable_model_postdose2_loglinear)

multivariable_model_postdose2_loglinear_allos = lm(formula = logic50~
                                                     #p values <0.1
                                                     age2+
                                                     COVID_predose2+
                                                     indication+
                                                     dosing_interval2+
                                                     time_since_HSCT2+
                                                     hastheparticipantreceivedrituxim+
                                                     no_prev_HSCT+
                                                     lymph_count_vac2+
                                                     #a priori
                                                     prime+
                                                     conditioning+#include now
                                                     lymph_depletion3,#include now
                                                   data = regression_WT_postdose2_allos)

summary(multivariable_model_postdose2_loglinear_allos)


####Healthy control linear regression#####

regression_WT_postdose2_HC <- subset(regression_WT_postdose2_HC,select = c(pid, ic50, prime, age, sex, dosing_interval, serostatus))

regression_WT_postdose2_HC$age2 <- ifelse(regression_WT_postdose2_HC$age < 45, "<45",
                                           ifelse(regression_WT_postdose2_HC$age < 65, "45-64", ">65"))

regression_WT_postdose2_HC$dosing_interval  <- ifelse(regression_WT_postdose2_HC$dosing_interval <42, "<6 weeks",
                                                     ifelse(regression_WT_postdose2_HC$dosing_interval >41, ">6 weeks","banana"))                

regression_WT_postdose2_HC$serostatus <- as.factor(regression_WT_postdose2_HC$serostatus)
regression_WT_postdose2_HC$serostatus <-relevel(regression_WT_postdose2_HC$serostatus, ref="N")

regression_WT_postdose2_HC$prime <- as.factor(regression_WT_postdose2_HC$prime)
regression_WT_postdose2_HC$prime <-relevel(regression_WT_postdose2_HC$prime, ref="mRNA")

regression_WT_postdose2_HC$logic50 <- log10(regression_WT_postdose2_HC$ic50)

summary(lm(logic50~prime, data=regression_WT_postdose2_HC))
summary(lm(logic50~age2, data=regression_WT_postdose2_HC))
summary(lm(logic50~sex, data=regression_WT_postdose2_HC))
summary(lm(logic50~serostatus, data=regression_WT_postdose2_HC))
summary(lm(logic50~dosing_interval, data=regression_WT_postdose2_HC))

lm_model_HC <- lm(formula = logic50~
                    age2+
                    sex+
                    serostatus+
                    prime+
                    dosing_interval,
                  data = regression_WT_postdose2_HC)

summary(lm_model_HC)

####Binding Antibody Regression####
###Roche Seroconversion regression post dose 2####

regression_S <- full_join(metadata,bindingab,by="pid")

regression_S <- regression_S %>%
  filter(serum_source=="HSCT")

regression_S$lymph_depletion3 <- ifelse(regression_S$lymph_depletion == "None", "No",
                                        ifelse(regression_S$lymph_depletion == "Anti-thymocyte globulin", "ATG",
                                               ifelse(regression_S$lymph_depletion == "Campath", "Non-ATG",
                                                      ifelse(regression_S$lymph_depletion == "Post-transplant Cyclophosphamide", "Non-ATG", 
                                                             ifelse(regression_S$lymph_depletion == "", "Unknown",
                                                                    ifelse(regression_S$lymph_depletion == "Unknown", "Unknown",
                                                                           regression_S$lymph_depletion))))))

regression_S$lymph_depletion3 <- as.factor(regression_S$lymph_depletion3)
regression_S$lymph_depletion3 <-relevel(regression_S$lymph_depletion3, ref="Non-ATG")

regression_S$age2 <- ifelse(regression_S$age < 45, "<45",
                            ifelse(regression_S$age < 65, "45-64", ">65"))

regression_S$visit_date <- as.Date(regression_S$visit_date)
regression_S$D0 <- as.Date(regression_S$D0)
regression_S$V1 <- as.Date(regression_S$V1)

regression_S$time_since_HSCT <- difftime(regression_S$V1, regression_S$D0, unit="days")

regression_S$time_since_HSCT2 <- ifelse(regression_S$time_since_HSCT <366, "<12 months",
                                        ifelse(regression_S$time_since_HSCT >365, ">12 months", "banana"))                

regression_S$time_since_HSCT2 <- as.factor(regression_S$time_since_HSCT2)
regression_S$time_since_HSCT2 <-relevel(regression_S$time_since_HSCT2, ref=">12 months")

regression_S$days_since_last_vac <- difftime(regression_S$visit_date, regression_S$V2, unit="days")

regression_S$ivig_date <- as.Date(regression_S$ivig_date, format="%d/%m/%Y")
regression_S$time_since_IVIG <- difftime(regression_S$visit_date, regression_S$ivig_date, unit="days")
regression_S <- regression_S %>%
  filter(time_since_IVIG > 90 | is.na(time_since_IVIG)) #gets rid of any participants who had IVIG within 90 days of postdose 2 visit

regression_S$dosing_interval2 <- ifelse(regression_S$dosing_interval <42, "<6 weeks",
                                        ifelse(regression_S$dosing_interval >41, ">6 weeks","banana"))                

regression_S$bindingab_vaccine_response <- as.factor(regression_S$bindingab_vaccine_response)

regression_S$bindingab_vaccine_response_normal <- ifelse(regression_S$bindingab_vaccine_response == "DE (normal)", "1", "0")

regression_S$bindingab_vaccine_response_detectable <- ifelse(regression_S$bindingab_vaccine_response == "DE (low)", "1", ###change this to 0 if wanting to group low/non responders together
                                                   ifelse(regression_S$bindingab_vaccine_response == "DE (normal)", "1",
                                                                 ifelse(regression_S$bindingab_vaccine_response == "NDE", "0", "banana")))

regression_S$bindingab_vaccine_response_normal <- as.factor(regression_S$bindingab_vaccine_response_normal)
regression_S$bindingab_vaccine_response_detectable <- as.factor(regression_S$bindingab_vaccine_response_detectable)

regression_S$prime <- as.factor(regression_S$prime)
regression_S$prime <-relevel(regression_S$prime, ref="mRNA")

regression_S$indication <- as.factor(regression_S$indication)
regression_S$indication <-relevel(regression_S$indication, ref="Myeloid")

regression_S$conditioning <- as.factor(regression_S$conditioning)
regression_S$conditioning <-relevel(regression_S$conditioning, ref="Reduced Intensity")

#univariate logistic regression for normal binding ab vaccine response####

lapply(c("prime", "age","age2", "sex",
         "time_since_HSCT2","days_since_last_vac",
         "dosing_interval2", "COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion3", "tbi",
         "no_prev_HSCT", "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"),
       function(var) {
         
         formula = as.formula(paste("bindingab_vaccine_response_normal~", var))
         
         res.logit = glm(formula, family=binomial(link="logit"),  data=regression_S)
         
         summary(res.logit)
         
       })


lapply(c("prime", "age","age2", "sex",
         "time_since_HSCT2","days_since_last_vac",
         "dosing_interval2", "COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion3", "tbi",
         "no_prev_HSCT", "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"),
       function(var) {
         
         formula = as.formula(paste("bindingab_vaccine_response_normal~", var))
         
         res.logit = glm(formula, family=binomial(link="logit"), data=regression_S)
         
         exp(coef(res.logit))
         
       })





#Get CIs for OR

lapply(c("prime", "age","age2", "sex",
         "time_since_HSCT2","days_since_last_vac",
         "dosing_interval2", "COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion3", "tbi",
         "no_prev_HSCT", "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "any_gvhd_vac1","any_gvhd_vac2"),
       
       function(var) {
         
         formula = as.formula(paste("bindingab_vaccine_response_normal~", var))
         
         res.logit = glm(formula, family=binomial(link="logit"), data=regression_S)
         
         exp(confint(res.logit))
         
       })


#univariate linear regression models for raw value and log value
summary(lm(PostBoostRocheS~prime, data=regression_S))
summary(lm(PostBoostRocheS~age2, data=regression_S))
summary(lm(PostBoostRocheS~sex, data=regression_S))
summary(lm(PostBoostRocheS~ethnicity, data=regression_S))
summary(lm(PostBoostRocheS~time_since_HSCT2, data=regression_S))
summary(lm(PostBoostRocheS~days_since_last_vac, data=regression_S))
summary(lm(PostBoostRocheS~HSCT_group, data=regression_S)) 
summary(lm(PostBoostRocheS~indication, data=regression_S))
summary(lm(PostBoostRocheS~diseasestatusattrialentry, data=regression_S))
summary(lm(PostBoostRocheS~conditioning, data=regression_S))
summary(lm(PostBoostRocheS~lymph_depletion, data=regression_S))
summary(lm(PostBoostRocheS~lymph_depletion3, data=regression_S))
summary(lm(PostBoostRocheS~tbi, data=regression_S))
summary(lm(PostBoostRocheS~prev_auto, data=regression_S))
summary(lm(PostBoostRocheS~prev_allo, data=regression_S))
summary(lm(PostBoostRocheS~no_prev_HSCT, data=regression_S))
#summary(lm(PostBoostRocheS~toci_rec, data=regression_S))
summary(lm(PostBoostRocheS~hastheparticipantreceivedrituxim, data=regression_S))
summary(lm(PostBoostRocheS~lymph_count_vac2, data=regression_S))
summary(lm(PostBoostRocheS~igg_level_vac2, data=regression_S))
summary(lm(PostBoostRocheS~lymph_infusion_vac2, data=regression_S))
summary(lm(PostBoostRocheS~hsct_hadnolymphstmcllpriorvac, data=regression_S))
summary(lm(PostBoostRocheS~COVID_predose2, data=regression_S))
summary(lm(PostBoostRocheS~any_gvhd_vac1, data=regression_S))
summary(lm(PostBoostRocheS~any_gvhd_vac2, data=regression_S))
summary(lm(PostBoostRocheS~dosing_interval2, data=regression_S))

#multivariate log linear model
regression_S_linear <- regression_S
regression_S_linear$PostBoostRocheS <- ifelse(regression_S_linear$PostBoostRocheS == 0, 1, regression_S_linear$PostBoostRocheS)
regression_S_linear$PostBoostRocheS_log10 <- log10(regression_S_linear$PostBoostRocheS)

multivariable_model_postdose2_Roche_linear = lm(formula = PostBoostRocheS_log10~
                                                  #p values <0.1
                                                  age2+
                                                  prime+
                                                  lymph_count_vac2+
                                                  time_since_HSCT2+
                                                  any_gvhd_vac1+
                                                  COVID_predose2+
                                                  #a priori
                                                  conditioning+
                                                  HSCT_group+ 
                                                  dosing_interval2+
                                                  indication+
                                                  hastheparticipantreceivedrituxim
                                                , data = regression_S_linear)

summary(multivariable_model_postdose2_Roche_linear)

####Multivariate logistic model for normal post-vacc binding ab response
multivariable_model_postdose2_Roche = glm(bindingab_vaccine_response_normal~
                                            #p values <0.1
                                            lymph_count_vac2+
                                            age2+
                                            dosing_interval2+
                                            time_since_HSCT2+
                                            hastheparticipantreceivedrituxim+
                                            #a priori
                                            HSCT_group+ 
                                            prime+ 
                                            indication+
                                            COVID_predose2,
                                          family=binomial(link="logit"), data = regression_S_linear)

#binding ab regression models for allos only
regression_S2_allos <- regression_S_linear %>%
  filter(HSCT_group=="Allo")

multivariable_model_postdose2_Roche_allos = glm(bindingab_vaccine_response_normal~
                                                  #p values <0.1
                                                  age2+
                                                  lymph_count_vac2+
                                                  time_since_HSCT2+
                                                  hastheparticipantreceivedrituxim+
                                                  #a priori
                                                  conditioning+ #only have data for allos
                                                  dosing_interval2+
                                                  prime+ 
                                                  indication+
                                                  COVID_predose2+
                                                  lymph_depletion3,#only have data for allos
                                                family=binomial(link="logit"), data = regression_S2_allos)

summary(multivariable_model_postdose2_Roche_allos)

###ELISpot regression####
elispot$panel_14 <- as.numeric(elispot$panel_14)
elispot$nil_control <- as.numeric(elispot$nil_control)
elispot$transformed_result <- elispot$panel_14 - elispot$nil_control 
    #Transformed Assay Results = (Raw Count - Nil Control) x 4 #units for transformed result are SFC/10^6 PBMC
    #Raw Count = Panel 14 result if participant is not recruited to the Renal Cohort
    #For this assay spot value of 4 or less is considered to be a non-response, the results provided from OI contain the assay results in continuous format and as such this definition will be used within the statistical analysis to categorise patients as responders/non-responders where relevant. 

elispot$final_result <- ifelse(elispot$transformed_result <= 4, "non-responder",
                               ifelse(elispot$transformed_result > 4, "responder", "banana"))

regression_elispot <- full_join(metadata,elispot,by="pid")

regression_elispot <- regression_elispot %>%
  filter(serum_source=="HSCT")

regression_elispot$ivig_date <- as.Date(regression_elispot$ivig_date)
regression_elispot$time_since_HSCT <- difftime(regression_elispot$V1, regression_elispot$D0, unit="days")

regression_elispot$time_since_HSCT2 <- ifelse(regression_elispot$time_since_HSCT <366, "<12 months",
                                               ifelse(regression_elispot$time_since_HSCT >365, ">12 months", "banana"))                

regression_elispot$time_since_HSCT2 <- as.factor(regression_elispot$time_since_HSCT2)
regression_elispot$time_since_HSCT2 <-relevel(regression_elispot$time_since_HSCT2, ref=">12 months")

regression_elispot$days_since_last_vac <- difftime(regression_elispot$visit_date, regression_elispot$V2, unit="days")
regression_elispot$igg_level_vac2 <- as.numeric(regression_elispot$igg_level_vac2)
regression_elispot$lymph_count_vac2 <- as.numeric(regression_elispot$lymph_count_vac2)

regression_elispot$final_result <- ifelse(regression_elispot$final_result == "responder", "1",
                                           ifelse(regression_elispot$final_result == "non-responder", "0", "banana"))

regression_elispot$final_result <- as.factor(regression_elispot$final_result)

regression_elispot$indication <- as.factor(regression_elispot$indication)
regression_elispot$indication <-relevel(regression_elispot$indication, ref="Myeloid")

regression_elispot$conditioning <- as.factor(regression_elispot$conditioning)
regression_elispot$conditioning <-relevel(regression_elispot$conditioning, ref="Reduced Intensity")

regression_elispot$hsct_hadnolymphstmcllpriorvac <- as.factor(regression_elispot$hsct_hadnolymphstmcllpriorvac)
regression_elispot$hsct_hadnolymphstmcllpriorvac <-relevel(regression_elispot$hsct_hadnolymphstmcllpriorvac, ref="Checked")

regression_elispot$hsct_hadstemcelltopuppriorvac <- as.factor(regression_elispot$hsct_hadstemcelltopuppriorvac)
regression_elispot$hsct_hadstemcelltopuppriorvac <-relevel(regression_elispot$hsct_hadstemcelltopuppriorvac, ref="Unchecked")

regression_elispot$hsct_hadlymphocyteinfuspriorvac <- as.factor(regression_elispot$hsct_hadlymphocyteinfuspriorvac)
regression_elispot$hsct_hadlymphocyteinfuspriorvac <- relevel(regression_elispot$hsct_hadlymphocyteinfuspriorvac, ref="Unchecked")

regression_elispot$prime <- as.factor(regression_elispot$prime)
regression_elispot$prime <-relevel(regression_elispot$prime, ref="mRNA")

regression_elispot$age2 <- ifelse(regression_elispot$age < 45, "<45",
                                   ifelse(regression_elispot$age < 65, "45-64", ">65"))

regression_elispot$dosing_interval2 <- ifelse(regression_elispot$dosing_interval <42, "<6 weeks",
                                               ifelse(regression_elispot$dosing_interval >41, ">6 weeks","banana"))   


lapply(c("prime", "age2", "sex",
         "time_since_HSCT2", "days_since_last_vac", "dosing_interval2","COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion", "tbi",
         "prev_auto", "prev_allo", "no_prev_HSCT", 
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "hsct_hadlymphocyteinfuspriorvac", "hsct_hadnolymphstmcllpriorvac",
         "any_gvhd_vac2"
),
function(var) {
  
  formula = as.formula(paste("final_result~", var))
  
  res.logit = glm(formula, family=binomial(link="logit"),  data=regression_elispot)
  
  summary(res.logit)
  
})


lapply(c("prime", "age2", "sex",
         "time_since_HSCT2", "days_since_last_vac", "dosing_interval2","COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion", "tbi",
         "prev_auto", "prev_allo", "no_prev_HSCT", 
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "hsct_hadlymphocyteinfuspriorvac", "hsct_hadnolymphstmcllpriorvac",
         "any_gvhd_vac2"),
       
       function(var) {
         
         formula = as.formula(paste("final_result~", var))
         
         res.logit = glm(formula, family=binomial(link="logit"), data=regression_elispot)
         
         exp(coef(res.logit))
         
       })





#Get CIs for OR

lapply(c("prime", "age2", "sex",
         "time_since_HSCT2", "days_since_last_vac", "dosing_interval2","COVID_predose2",
         "HSCT_group", "indication", "diseasestatusattrialentry",
         "conditioning","lymph_depletion", "tbi",
         "prev_auto", "prev_allo", "no_prev_HSCT", 
         "hastheparticipantreceivedrituxim", 
         "lymph_count_vac2", "igg_level_vac2", 
         "hsct_hadlymphocyteinfuspriorvac", "hsct_hadnolymphstmcllpriorvac",
         "any_gvhd_vac2"),
       
       function(var) {
         
         formula = as.formula(paste("final_result~", var))
         
         res.logit = glm(formula, family=binomial(link="logit"), data=regression_elispot)
         
         exp(confint(res.logit))
         
       })

