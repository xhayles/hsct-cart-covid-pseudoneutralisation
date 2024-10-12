
regression_elispot$transformed_result <- as.numeric(regression_elispot$transformed_result)

regression_elispot3 <- regression_elispot %>%
  filter(prime == "mRNA" | prime == "AZ")

wilcox_test(transformed_result ~ prime, data=regression_elispot3) %>%
  adjust_pvalue(method = "bonferroni") 

figure3_supp <- ggplot(regression_elispot3,aes(x=prime, y=transformed_result, color=prime, fill=prime))+
  geom_jitter()+
  geom_boxplot(outlier.shape=NA, alpha=0.1)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(expression(paste("SFC/10"^"6"*" PBMC")))+
  xlab("Vaccination Received for V1+V2")+
  geom_hline(yintercept = 4, linetype="dashed", color="grey")+
  scale_y_break(c(130,480))+
  scale_y_continuous(limits=c(0,500), breaks=c(0,25,50,75,100,125,480,500))+
  scale_color_manual(values=c( "#007d79","#6929c4"))+
  scale_fill_manual(values=c( "#007d79","#6929c4"))+
  annotate("text", label="p=0.584",x = 1.5, y = 130,size=4, color="black")
