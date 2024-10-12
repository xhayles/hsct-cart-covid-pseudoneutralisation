
all_ic50_results_v2 <- all_ic50_results3 %>%
  filter(visit_code == "Post-V2")

all_ic50_results_v2$logic50 <- log10(all_ic50_results_v2$ic50)

all_ic50_results_v2_N <- all_ic50_results_v2 %>%
  filter(serostatus == "N")

ic50_naive_medians <- all_ic50_results_v2_N %>%
  group_by(pseudovirus, prime,serum_source)%>%
  dplyr::summarise(median=median(ic50))

ic50_naive_medians <- ic50_naive_medians %>%
  pivot_wider(names_from = 'serum_source', values_from = 'median')

fold_change_naive<- ic50_naive_medians %>%
  mutate(foldChange = HC/HSCT)

ann_text_naive<- fold_change_naive
ann_text_naive$serum_source = "HSCT"
ann_text_naive$visit_code = "Post-V2"
ann_text_naive$ic50 = 500000
ann_text_naive$logic50 = 6.5
ann_text_naive$foldChange <- round(ann_text_naive$foldChange, digits = 2)


stat.testN <- all_ic50_results_v2_N %>%
  group_by(pseudovirus,prime) %>%
  wilcox_test(logic50 ~ serum_source) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") #calculate p values between the groups

stat.testN <- stat.testN %>%
  add_y_position() %>%
  add_x_position()

stat.testN$xmin <- 0.8 #shifts the brackets and p values to the correct location
stat.testN$xmax <- 1.2 #shifts the brackets and p values to the correct location

stat.testN$p.adj <- ifelse(stat.testN$p.adj < 0.001, "p <0.001", paste0("p = ", round(stat.testN$p.adj, 3)))

figure4c <- ggplot(all_ic50_results_v2_N, aes(x=visit_code, y=log10(ic50), color=serum_source))+
  facet_nested(prime~pseudovirus,
               labeller = labeller(pseudovirus = pv.labs))+
  geom_boxplot(outlier.shape=NA, coef=Inf) + 
  theme_bw()+
  geom_point(position=position_jitterdodge(), alpha=0.3, size=1)+
  scale_color_manual(values=c("#000000","#0041C2"), labels=c("Healthy Controls", "HSCT/CAR-T"))+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_text(face="bold"))+
  xlab(NULL)+
  ylab("log10 NT50")+
  scale_y_continuous(limits = c(0,7), breaks=c(0,1,2,3,4,5))+
  ggtitle("C)")+
  geom_hline(yintercept=1.60206, linetype="dashed", color="black")+
  geom_hline(yintercept=6, linetype="solid", color="grey")+
  stat_pvalue_manual(stat.testN, label="p.adj", size=3)+
  geom_text(data = ann_text_naive, aes(label = paste(foldChange, "x")), x = 1, y = 6.5,size=3, color="black")

#supplementary figure 2
all_ic50_results_v2_C <- all_ic50_results_v2 %>%
  filter(serostatus == "C")

ic50_convalescent_medians <- all_ic50_results_v2_C %>%
  group_by(pseudovirus, prime,serum_source)%>%
  dplyr::summarise(median=median(ic50))

ic50_convalescent_medians <- ic50_convalescent_medians %>%
  pivot_wider(names_from = 'serum_source', values_from = 'median')

fold_change_conv<- ic50_convalescent_medians %>%
  mutate(foldChange = HC/HSCT)

stat.testC <- all_ic50_results_v2_C %>%
  group_by(pseudovirus,prime) %>%
  wilcox_test(logic50 ~ serum_source) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") #calculate p values between the groups

stat.testC <- stat.testC %>%
  add_y_position() %>%
  add_x_position()

stat.testC$xmin <- 0.8 #shifts the brackets and p values to the correct location
stat.testC$xmax <- 1.2 #shifts the brackets and p values to the correct location

stat.testC$p.adj <- ifelse(stat.testC$p.adj < 0.001, "p <0.001", paste0("p = ", round(stat.testC$p.adj, 3)))

ann_text_conv <- fold_change_conv
ann_text_conv$serum_source = "HSCT"
ann_text_conv$visit_code = "Post-V2"
ann_text_conv$ic50 = 500000
ann_text_conv$logic50 = 6.5
ann_text_conv$foldChange <- round(ann_text_conv$foldChange, digits = 2)

figure2_supp <- ggplot(all_ic50_results_v2_C, aes(x=visit_code, y=log10(ic50), color=serum_source))+
  facet_grid(prime~pseudovirus,
             labeller = labeller(pseudovirus = pv.labs))+
  geom_boxplot(outlier.shape=NA, coef=Inf)+
  theme_bw()+
  geom_point(position=position_jitterdodge(), alpha=0.3, size=1)+
  scale_color_manual(values=c("#000000","#0041C2"), labels=c("Healthy Controls", "HSCT/CAR-T"))+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab(NULL)+
  ylab("log10 NT50")+
  scale_y_continuous(limits = c(0,5.7), breaks=c(0,1,2,3,4,5))+
  ggtitle(NULL)+
  geom_hline(yintercept=1.60206, linetype="dashed", color="black")+
  geom_hline(yintercept=5.5, linetype="solid", color="grey")+
  stat_pvalue_manual(stat.testC, label="p.adj", size=3)+
  geom_text(data = ann_text_conv, aes(label = paste(foldChange, "x")), x = 1, y = 5.7,size=3, color="black")


#Figure 4a
rcdf <- subset(all_ic50_results4, select=c(sample.id,ic50,pseudovirus,visit_code,serum_source,prime,serostatus))
rcdf$group <- paste(rcdf$serostatus,"_",rcdf$prime)
rcdf_wt <- rcdf %>%
  filter(pseudovirus=="WT" & visit_code == "Post-V2")

rcdf_wt <- subset(rcdf_wt, select=c(ic50,group,serum_source))

mrcdf_wt <- melt(rcdf_wt, id = c("group","serum_source"))

mrcdf_wt$group2 <- paste(mrcdf_wt$serum_source,mrcdf_wt$group)

mrcdf_wt$group <- ifelse(mrcdf_wt$group == "N _ AZ", "AZ (Naive)",
                         ifelse(mrcdf_wt$group == "C _ AZ", "AZ (Convalescent)",
                                ifelse(mrcdf_wt$group == "N _ mRNA", "mRNA (Naive)",
                                       ifelse(mrcdf_wt$group == "C _ mRNA", "mRNA (Convalescent)", "banana"))))

mrcdf_wt$group <- factor(mrcdf_wt$group, levels=c("AZ (Naive)",  "AZ (Convalescent)", "mRNA (Naive)", "mRNA (Convalescent)"))

serum_source.labs <- c("Healthy Controls", "HSCT/CAR-T")
names(serum_source.labs) <- c("HC", "HSCT")

figure4a <- ggplot(mrcdf_wt, aes(log10(value), colour = group)) +
  facet_grid(serum_source~.,
             labeller = labeller(serum_source = serum_source.labs))+
  stat_ecdf(geom = "step", linewidth=1) +
  theme_bw()+
  theme(legend.position = "top",
        plot.title=element_text(face="bold"))+
  scale_x_continuous(name="log10 NT50") +
  scale_y_reverse(name="Reverse Cumulative Distribution (%)", breaks=c(1.00, 0.75, 0.50, 0.25, 0.00), labels=c("0", "25", "50", "75", "100")) +
  ggtitle(NULL)+
  scale_color_manual(values=c("#be95ff", #N_AZ
                              "#6929c4",#C_AZ
                              "#3ddbd9", #N_mRNA,
                              "#007d79"#C_mRNA
  ))+
  geom_vline(xintercept=log10(40), linetype="dashed")+
  ggtitle("A)")
