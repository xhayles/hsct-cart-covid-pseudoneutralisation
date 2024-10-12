
HC <- exp(coef(lm_model_HC)["primeAZ"])
HSCT <- exp(coef(multivariable_model_postdose2_loglinear)["primeAZ"])
HC_95CI <- exp(confint(lm_model_HC)["primeAZ", ])
HSCT_95CI <- exp(confint(multivariable_model_postdose2_loglinear)["primeAZ", ])

iteration_result2 <- data.frame(group=c("HC","HSCT"),
                                adj_diff_mean = c(HC, HSCT),
                                lower_CI = c(HC_95CI[1], HSCT_95CI[1]),
                                upper_CI = c(HC_95CI[2],HSCT_95CI[2]))

iteration_result2$group <- factor(iteration_result2$group,
                                  levels=c("HSCT",
                                           "HC"))

serum_source.labs <- c("Healthy Controls", "HSCT/CAR-T")
names(serum_source.labs) <- c("HC", "HSCT")

figure4b<- ggplot(iteration_result2, aes(x=adj_diff_mean, y=group))+
  geom_point(color="black")+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbar(aes(xmin=lower_CI, xmax=upper_CI), width=0, color="black")+
  theme_bw()+
  theme(axis.title.x = element_text(hjust = 0.3),
        plot.title=element_text(face="bold"))+
  xlab("\n
       Geometric Mean Ratio of log10 NT50")+
  ylab(NULL)+
  scale_x_continuous(breaks=c(0.5,0.75,1,1.25,1.5,1.75,2.00))+
  scale_y_discrete(position="left", labels=serum_source.labs)+
  ggtitle("B)")

figure4b_arrow <-
  figure4b + annotate("segment", x =   1.1 , xend = 2.05, y = 0, yend = 0,colour = "black", arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = 1.2, y = -0.2, label = "favours AZ", size=3.5) +
  coord_cartesian(xlim = c(0.5, 2),  ylim=c(1, 2), clip="off")+
  annotate("segment", x =   0.9 , xend = 0.45, y = 0 , yend = 0, colour = "black", arrow = arrow(length = unit(0.2, "cm")))+
  annotate("text", x = 0.8, y = -0.2, label = "favours mRNA", size=3.5)
