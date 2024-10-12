#to create dataframe for Figure 2b
`factor_<45` <- 1
`factor_45-64` <- exp(coef(multivariable_model_postdose2_loglinear)["age245-64"])
`CI_45-64`<- exp(confint(multivariable_model_postdose2_loglinear)["age245-64", ])
`factor_>65` <- exp(coef(multivariable_model_postdose2_loglinear)["age2>65"])
`CI_>65`<- exp(confint(multivariable_model_postdose2_loglinear)["age2>65", ])

`factor_covid_no` <- 1
`factor_covid_yes` <- exp(coef(multivariable_model_postdose2_loglinear)["COVID_predose2Yes"])
`CI_covid_yes`<- exp(confint(multivariable_model_postdose2_loglinear)["COVID_predose2Yes", ])

#`factor_RIC` <- 1
#`factor_myeloablative` <- exp(coef(multivariable_model_postdose2_loglinear_allos)["conditioningMyeloablative"])
#`CI_myeloablative`<- exp(confint(multivariable_model_postdose2_loglinear_allos)["conditioningMyeloablative", ])
#`factor_nonmyeloablative` <- exp(coef(multivariable_model_postdose2_loglinear_allos)["conditioningNon-myeloablative"])
#`CI_nonmyeloablative`<- exp(confint(multivariable_model_postdose2_loglinear_allos)["conditioningNon-myeloablative", ])

`factor_noprevHSCT` <- 1
`factor_prevHSCT` <- exp(coef(multivariable_model_postdose2_loglinear)["no_prev_HSCTUnchecked"])
`CI_prevHSCT`<- exp(confint(multivariable_model_postdose2_loglinear)["no_prev_HSCTUnchecked", ])

`factor_ritux_no` <- 1
`factor_ritux_yes` <- exp(coef(multivariable_model_postdose2_loglinear)["hastheparticipantreceivedrituximYes"])
`CI_ritux_yes`<- exp(confint(multivariable_model_postdose2_loglinear)["hastheparticipantreceivedrituximYes", ])

`factor_prime_mRNA` <- 1
`factor_prime_AZ` <- exp(coef(multivariable_model_postdose2_loglinear)["primeAZ"])
`CI_prime_AZ`<- exp(confint(multivariable_model_postdose2_loglinear)["primeAZ", ])
#`factor_prime_mixed` <- exp(coef(multivariable_model_postdose2_loglinear)["primeMixed"])
#`CI_prime_mixed`<- exp(confint(multivariable_model_postdose2_loglinear)["primeMixed", ])

`factor_HSCT>12m` <- 1
`factor_HSCT<12m` <- exp(coef(multivariable_model_postdose2_loglinear)["time_since_HSCT2<12 months"])
`CI_HSCT<12m`<- exp(confint(multivariable_model_postdose2_loglinear)["time_since_HSCT2<12 months", ])

`factor_dosing<6` <- 1
`factor_dosing>6` <- exp(coef(multivariable_model_postdose2_loglinear)["dosing_interval2>6 weeks"])
`CI_dosing>6`<- exp(confint(multivariable_model_postdose2_loglinear)["dosing_interval2>6 weeks", ])

`factor_myeloid` <- 1
`factor_lymphoid` <- exp(coef(multivariable_model_postdose2_loglinear)["indicationLymphoid"])
`CI_lymphoid`<- exp(confint(multivariable_model_postdose2_loglinear)["indicationLymphoid", ])

`factor_auto` <- 1
`factor_allo` <- exp(coef(multivariable_model_postdose2_loglinear)["HSCT_groupAllo"])
`CI_allo`<- exp(confint(multivariable_model_postdose2_loglinear)["HSCT_groupAllo", ])
#`factor_cart` <- exp(coef(multivariable_model_postdose2_loglinear)["HSCT_groupCAR-T"])
#`CI_cart`<- exp(confint(multivariable_model_postdose2_loglinear)["HSCT_groupCAR-T", ])

#`factor_nonatg` <- 1
#`factor_atg` <- exp(coef(multivariable_model_postdose2_loglinear_allos)["lymph_depletion3ATG"])
#`CI_atg`<- exp(confint(multivariable_model_postdose2_loglinear_allos)["lymph_depletion3ATG", ])
#`factor_none` <- exp(coef(multivariable_model_postdose2_loglinear_allos)["lymph_depletion3No"])
#`CI_none`<- exp(confint(multivariable_model_postdose2_loglinear_allos)["lymph_depletion3No", ])

iteration_result3 <- data.frame(group=c("Age","Age","Age",
                                        "Previous COVID","Previous COVID",
                                        "Previous HSCT", "Previous HSCT",
                                        "Rituximab","Rituximab",
                                        "Vaccine", "Vaccine", 
                                        "HSCT D0", "HSCT D0",
                                        "Indication", "Indication",
                                        "HSCT", "HSCT", 
                                        "Dosing interval", "Dosing interval"),
                                group2=c("<45","45-64","65+",
                                         "No","Yes",
                                         "No", "Yes",
                                         "No","Yes",
                                         "mRNA", "AZ",
                                         ">12 months", "<12 months",
                                         "Myeloid", "Lymphoid",
                                         "Auto","Allo",
                                         "<6 weeks", ">6 weeks"),
                                adj_diff_mean=c(`factor_<45`,`factor_45-64`,`factor_>65`,
                                                `factor_covid_no` ,`factor_covid_yes`,
                                                `factor_noprevHSCT`,`factor_prevHSCT`,
                                                `factor_ritux_no`,`factor_ritux_yes`,
                                                `factor_prime_mRNA`,`factor_prime_AZ`,
                                                `factor_HSCT>12m`,`factor_HSCT<12m`,
                                                `factor_myeloid`,`factor_lymphoid`,
                                                `factor_auto`,`factor_allo`,
                                                `factor_dosing<6`,`factor_dosing>6`),
                                lower_CI =c(1,`CI_45-64`[1],`CI_>65`[1],
                                            1,`CI_covid_yes`[1],
                                            1,`CI_prevHSCT`[1],
                                            1,`CI_ritux_yes`[1],
                                            1,`CI_prime_AZ`[1],
                                            1,`CI_HSCT<12m`[1],
                                            1,`CI_lymphoid`[1],
                                            1,`CI_allo`[1],
                                            1, `CI_dosing>6`[1]),
                                upper_CI =c(1,`CI_45-64`[2],`CI_>65`[2],
                                            1,`CI_covid_yes`[2],
                                            1,`CI_prevHSCT`[2],
                                            1,`CI_ritux_yes`[2],
                                            1,`CI_prime_AZ`[2],
                                            1,`CI_HSCT<12m`[2],
                                            1,`CI_lymphoid`[2],
                                            1,`CI_allo`[2],
                                            1, `CI_dosing>6`[2]),
                                signify.n=c(0,0,0, #manually or 1 for significance, will change colour of line/point
                                            0,0,
                                            0,0,
                                            0,0,
                                            0,0,
                                            0,0,
                                            0,0,
                                            0,0,
                                            0,0))

iteration_result3$group3 <- factor(iteration_result3$group2,
                                   levels=c("65+","45-64","<45",
                                            "Yes","No",
                                            "<12 months",">12 months",
                                            "Allo", "Auto",
                                            "Lymphoid", "Myeloid",
                                            "AZ","mRNA",
                                            ">6 weeks", "<6 weeks"))

iteration_result3$group <- factor(iteration_result3$group,
                                  levels=c(
                                    "Age",
                                    "Previous COVID",
                                    "HSCT",
                                    "HSCT D0",
                                    "Indication", 
                                    "Rituximab",
                                    "Vaccine",
                                    "Dosing interval",
                                    "Previous HSCT"
                                  ))

#data input file if preferred
iteration_result3 <- read.csv("./data/linear_regression_output.csv") 

figure3b<-
  ggplot(iteration_result3, aes(x=adj_diff_mean, y=group3,colour = signify.n))+
  facet_wrap(~group,
             ncol=1,
             scales="free_y",
             labeller = as_labeller(c(Rituximab="Rituximab",
                                      `Previous HSCT` ="Prev HSCT",
                                      `Previous COVID`="Prior COVID-19",
                                      #Conditioning="Conditioning*",
                                      Age ="Age (years)",
                                      Vaccine ="Vaccine", 
                                      `HSCT D0`="D0",
                                      Indication = "Indication", 
                                      HSCT = "HSCT/CAR-T",
                                      # `Lymph depletion` = "Lymph depletion*",
                                      `Dosing interval` = "Interval")),
             strip.position="left")+
  geom_point()+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbar(aes(xmin=lower_CI, xmax=upper_CI), width=0)+
  theme_classic()+
  xlab("Geometric Mean Ratio of log10 NT50")+
  ylab(NULL)+
  scale_y_discrete(position="left")+
  theme(legend.position = "none",
        axis.text.y=element_text(size=8),
        strip.text.y = element_text(size=6.5, face="bold"),
        strip.background = element_rect(colour="black", fill=NA),
        strip.placement = "outside",
        plot.title = element_text(hjust = -0.172))+
  scale_x_log10(breaks=c(0.2,0.5,1,2,5),labels=label_number(drop0trailing=TRUE))+
  ggtitle(" ")#space to allow ggarrange labels

###figure3a####
#data input file
OCTAVE_lownovshigh_model <- read.csv("./data/log_regression_output.csv")

OCTAVE_lownovshigh_model$Variable <- factor(OCTAVE_lownovshigh_model$Variable,
                                            levels=c("65+","45-64","<45",
                                                     "Yes","No",
                                                     "CAR-T", "Allo","Auto",
                                                     "<12 months",">12 months",
                                                     "Lymphoid", "Myeloid",
                                                     "Mixed","AZ","mRNA",
                                                     ">6 weeks", "<6 weeks",
                                                     "Count"
                                            ))

OCTAVE_lownovshigh_model$group <- factor(OCTAVE_lownovshigh_model$group,
                                         levels=c("Age",
                                                  "COVID",
                                                  "HSCT",
                                                  "HSCT D0",
                                                  "Indication",
                                                  "Rituximab",
                                                  "Vaccine",
                                                  "Dosing interval",
                                                  "Prev HSCT",
                                                  "Lymphocyte"))

figure3a <- 
  ggplot(OCTAVE_lownovshigh_model, aes(x=OR, y=Variable,colour = signify.n))+
  facet_wrap(~group,
             ncol=1,
             scales="free_y",
             labeller = as_labeller(c(Rituximab="Rituximab",
                                      Vaccine ="Vaccine",
                                      `COVID`="Prior COVID-19",
                                      #Conditioning="Conditioning*",
                                      #`Lymph depletion` = "Lymph depletion*"
                                      Age ="Age (years)",
                                      Indication = "Indication",
                                      `HSCT D0` = "D0",
                                      Lymphocyte = "Lymphocytes",
                                      `Dosing interval` = "Interval",
                                      HSCT = "HSCT/CAR-T",
                                      `Prev HSCT`="Prev HSCT")),
             strip.position="left")+
  geom_point()+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbar(aes(xmin=lower_CI, xmax=upper_CI), width=0)+
  theme_classic()+
  xlab("Odds Ratio of Neutralisation")+
  ylab(NULL)+
  scale_y_discrete(position="left")+
  theme(legend.position = "none",
        axis.text.y=element_text(size=8),
        strip.text.y = element_text(size=6.5, face="bold"),
        strip.background = element_rect(colour="black", fill=NA),
        strip.placement = "outside",
        plot.title = element_text(hjust = -0.172))+
  scale_color_manual(values=c("#13273b","#800020","#55b1f5"))+
  ggtitle(" ")+
  scale_x_log10(labels=label_number(drop0trailing=TRUE), breaks=log_breaks(7))
