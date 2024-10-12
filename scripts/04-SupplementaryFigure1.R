
regression_S$PostBoostRocheS_log <- ifelse(regression_S$PostBoostRocheS == 0, 0,
                                           log10(regression_S$PostBoostRocheS))

regression_S_dosingintervals <- subset(regression_S, select=c(pid,prime,dosing_interval,PostBoostRocheS_log, serostatus))
regression_S_dosingintervals <- regression_S_dosingintervals %>%
  filter(prime=="mRNA" | prime =="AZ")

labels <- regression_S_dosingintervals %>%
  dplyr::group_by(prime) %>%
  dplyr::summarise(r2 = cor(dosing_interval, PostBoostRocheS_log)^2)

labels$r2 <- sprintf("italic(R^2) == %.2f", labels$r2)

figure1_supp <- regression_S_dosingintervals %>%
  filter(prime == "mRNA" | prime == "AZ") %>%
  ggplot(aes(x=dosing_interval, y=PostBoostRocheS_log, color=prime))+
  facet_grid(prime~.)+
  geom_jitter(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(20,40,60,80,120,160,200,240))+ #removes 2 outliers
  scale_y_continuous()+
  scale_color_manual(values=c( "#007d79","#6929c4"))+
  xlab("Dosing interval, days")+
  ylab("RBD IgG at post-V2 visit, log10 AU/ml")+
  geom_hline(yintercept = log10(400), linetype="dashed", color="grey")+stat_cor(method = "pearson", label.x=100, label.y=2.9)
