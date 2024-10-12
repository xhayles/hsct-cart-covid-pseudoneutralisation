
all_ic50_results3<- regression

all_ic50_results3$serum_source <- ifelse(all_ic50_results3$serum_source =="PITCH", "HC",
                                         ifelse(all_ic50_results3$serum_source =="HERO", "NC",
                                                ifelse(all_ic50_results3$serum_source =="HSCT" & all_ic50_results3$visit_code == "postdose3", "HSCT (LR)",
                                                       ifelse(all_ic50_results3$serum_source =="HSCT (PROSECO)", "HSCT (NR)",
                                                              all_ic50_results3$serum_source))))

all_ic50_results3$serum_source2 <- ifelse(all_ic50_results3$visit_code =="neg_control", "Naive HC",
                                          ifelse(all_ic50_results3$visit_code =="baseline", "Convalescent HC", all_ic50_results3$serum_source))

all_ic50_results3$visit_code <- ifelse(all_ic50_results3$visit_code =="neg_control", "Vaccine-naive",
                                       ifelse(all_ic50_results3$visit_code =="baseline", "Vaccine-naive", 
                                              ifelse(all_ic50_results3$visit_code =="postdose1", "Post-V1", 
                                                     ifelse(all_ic50_results3$visit_code =="postdose2", "Post-V2", 
                                                            ifelse(all_ic50_results3$visit_code =="postdose3", "Post-V3", "banana")))))

all_ic50_results3$pseudovirus <- factor(all_ic50_results3$pseudovirus, levels=c("WT", "BA1", "BA5", "BQ", "XBB"))
all_ic50_results3$visit_code <- factor(all_ic50_results3$visit_code, levels=c("Vaccine-naive", "Post-V1", "Post-V2", "Post-V3"))

baselines <- all_ic50_results3 %>%
  filter(visit_code == "Vaccine-naive")

pv.labs <- c("Ancestral/B.1", "BA.1", "BA.5", "BQ.1.1", "XBB")
names(pv.labs) <- c("WT", "BA1", "BA5", "BQ", "XBB")

baselines$serum_source2 <- ifelse(baselines$serum_source2 == "Convalescent HC", "HC (C)",
                                  ifelse(baselines$serum_source2 == "Naive HC", "HC (N)",baselines$serum_source2))

figure1a <- ggplot(baselines, aes(x=visit_code, y=log10(ic50), color=serum_source2))+
  facet_nested(.~pseudovirus,
               labeller = labeller(pseudovirus = pv.labs))+
  geom_boxplot(outlier.shape=NA) + 
  theme_bw()+
  geom_point(position=position_jitterdodge(), alpha=0.3, size=1)+
  scale_color_manual(values=c("red","#636363"),
                     labels=c("Healthy Controls (convalescent)", "Healthy Controls (naive)"))+
  theme(legend.position = "right",
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_text(face="bold"))+
  xlab(NULL)+
  ylab("log10 NT50")+
  scale_y_continuous(limits = c(0,5))+
  ggtitle("A)")+
  geom_hline(yintercept=1.60206, linetype="dashed", color="black")

all_ic50_results4 <- all_ic50_results3 %>%
  filter(visit_code != "Vaccine-naive")

all_ic50_results4$logic50 <- log10(all_ic50_results4$ic50)

ic50_medians <- all_ic50_results4 %>%
  group_by(visit_code,pseudovirus, prime,serum_source)%>%
  dplyr::summarise(median=median(ic50))

ic50_medians <- ic50_medians %>%
  pivot_wider(names_from = 'serum_source', values_from = 'median')

ic50_medians_notV3 <- ic50_medians %>%
  filter(visit_code!= "Post-V3")

fc_median_notV3<- ic50_medians_notV3 %>%
  mutate(foldChange = HC/HSCT)

ic50_median_V3 <- ic50_medians %>%
  filter(visit_code== "Post-V3")

fc_median_V3<- ic50_median_V3 %>%
  mutate(foldChange = HC/`HSCT (LR)`) %>%
  mutate(foldChange2 = HC/`HSCT (NR)`) 

fc_median_notV3$foldChange2 <- NA

fc_medians<-rbind(fc_median_notV3, fc_median_V3)  #gives fold change SPLIT by prime for graph

###

ic50_medians_notsplit <- all_ic50_results4 %>%
  group_by(visit_code,pseudovirus,serum_source)%>%
  dplyr::summarise(median=median(ic50))

all_ic50_results4 %>%
  group_by(visit_code,pseudovirus,serum_source)%>%
  dplyr::summarise(min=min(ic50),max=max(ic50))

ic50_medians_notsplit <- ic50_medians_notsplit %>%
  pivot_wider(names_from = 'serum_source', values_from = 'median')

ic50_medians_notV3_notsplit <- ic50_medians_notsplit %>%
  filter(visit_code!= "Post-V3")

fc_median_notV3_notsplit<- ic50_medians_notV3_notsplit %>%
  mutate(foldChange = HC/HSCT)

ic50_median_V3_notsplit <- ic50_medians_notsplit %>%
  filter(visit_code== "Post-V3")

fc_median_V3_notsplit<- ic50_median_V3_notsplit %>%
  mutate(foldChange = HC/`HSCT (LR)`) %>%
  mutate(foldChange2 = HC/`HSCT (NR)`) 

fc_median_notV3_notsplit$foldChange2 <- NA

supp_table5_b<-rbind(fc_median_notV3_notsplit, fc_median_V3_notsplit) #gives fold change NOT split by prime

###

ann_text<- fc_medians #whichever fold change using here to annotate 
ann_text$ic50 = 500000
ann_text$logic50 = 6.5
ann_text$foldChange <- round(ann_text$foldChange, digits = 2)

ann_text2<- fc_medians %>% #whichever fold change using here to annotate 
  filter(!is.na(foldChange2))

ann_text2$logic50 = 7
ann_text2$foldChange2 <- round(ann_text2$foldChange2, digits = 2)

stat.test <- all_ic50_results4 %>%
  group_by(visit_code, pseudovirus,prime) %>%
  wilcox_test(logic50 ~ serum_source) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") #calculate p values between the groups

stat.test <- stat.test %>%
  add_y_position() %>%
  add_x_position(x="visit_code")

stat.test$x <- stat.test$x - 1 #shifts the brackets and p values to the correct location
stat.test$xmin <- stat.test$xmin - 1 #shifts the brackets and p values to the correct location
stat.test$xmax <- stat.test$xmax - 1 #shifts the brackets and p values to the correct location

stat.test$p.adj <- ifelse(stat.test$p.adj < 0.001, "p <0.001", paste0("p = ", round(stat.test$p.adj, 3)))

#splitting this into V1+2 and then V3

all_ic50_results4_a <- all_ic50_results4 %>%
  filter(visit_code != "Post-V3")

stat.test_a <- all_ic50_results4_a %>%
  group_by(visit_code, pseudovirus,prime) %>%
  wilcox_test(logic50 ~ serum_source) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") #calculate p values between the groups

stat.test_a <- stat.test_a %>%
  add_y_position() %>%
  add_x_position(x="visit_code")

stat.test_a$x <- stat.test_a$x - 1 #shifts the brackets and p values to the correct location
stat.test_a$xmin <- stat.test_a$xmin - 1 #shifts the brackets and p values to the correct location
stat.test_a$xmax <- stat.test_a$xmax - 1 #shifts the brackets and p values to the correct location

stat.test_a$p.adj <- ifelse(stat.test_a$p.adj < 0.001, "p <0.001", paste0("p = ", round(stat.test_a$p.adj, 3)))

ann_text_a <- ann_text %>%
  filter(visit_code != "Post-V3")

figure1b <- ggplot(all_ic50_results4_a, aes(x=visit_code, y=logic50, color=serum_source))+
  facet_nested(prime~pseudovirus,
               labeller = labeller(pseudovirus = pv.labs))+
  geom_boxplot(outlier.shape=NA, coef=Inf) + 
  theme_bw()+
  geom_point(position=position_jitterdodge(), alpha=0.3, size=1)+
  scale_color_manual(values=c("#000000","#0041C2"), labels=c("Healthy Controls", "HSCT/CAR-T"))+
  theme(legend.position = "right",
        legend.title=element_blank(),
        plot.title=element_text(face="bold"))+
  scale_y_continuous(limits = c(0,7), breaks=c(0,1,2,3,4,5))+
  xlab(NULL)+
  ylab("log10 NT50")+
  ggtitle("B)")+
  geom_hline(yintercept=1.60206, linetype="dashed", color="black")+
  geom_hline(yintercept=6, linetype="solid", color="grey")+
  stat_pvalue_manual(stat.test_a, label="p.adj", size=3)+
  geom_text(data = ann_text_a, aes(label =  paste(foldChange, "x")), size=3, color="black")

all_ic50_results4_b <- all_ic50_results4 %>%
  filter(visit_code == "Post-V3")


all_ic50_results4_b$prime <- ifelse(all_ic50_results4_b$prime == "AZ", "Prior AZ",
                                    ifelse(all_ic50_results4_b$prime == "mRNA", "Prior mRNA", all_ic50_results4_b$prime))

stat.test_b <- all_ic50_results4_b %>%
  group_by(visit_code, pseudovirus,prime) %>%
  wilcox_test(logic50 ~ serum_source) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") #calculate p values between the groups

stat.test_b <- stat.test_b %>%
  add_y_position() %>%
  add_x_position("prime")

stat.test_b$p.adj <- ifelse(stat.test_b$p.adj < 0.001, "p <0.001", paste0("p = ", round(stat.test_b$p.adj, 3)))

#stat.test_b <- stat.test_b %>%
# filter(group1 != "HSCT (LR)") #removes the comparison between LR and NR for post dose 3

ann_text_b <- ann_text %>%
  filter(visit_code == "Post-V3")
ann_text_b$logic50 = 7

ann_text_b$prime <- ifelse(ann_text_b$prime == "AZ", "Prior AZ",
                           ifelse(ann_text_b$prime == "mRNA", "Prior mRNA", ann_text_b$prime))

ann_text2_b <- ann_text2 %>%
  filter(visit_code == "Post-V3")


ann_text2_b$prime <- ifelse(ann_text2_b$prime == "AZ", "Prior AZ",
                            ifelse(ann_text2_b$prime == "mRNA", "Prior mRNA", ann_text2_b$prime))
ann_text2_b$logic50 = 6.5

all_ic50_results4_b$serum_source <- factor(all_ic50_results4_b$serum_source, levels=c("HC",  "HSCT (NR)", "HSCT (LR)"))


# Load the tibble package
library(tibble)

# Create the tibble manually based on fold changes
figure2_table <- data.frame(
  Group = c("","HSCT/CAR-T (normal responder)", "HSCT/CAR-T (low responder)"),
  `B.1/Ancestral` = c("Prior AZ", "0.73x", "4.62x"),
  ` ` = c("Prior mRNA","2.35x", "50.38x"),
  `BA.1` = c("Prior AZ","1.71x", "19.12x"),
  ` ` = c("Prior mRNA","4.59x", "16.44x"),
  `BA.5` = c("Prior AZ","0.47x", "7.09x"),
  ` ` = c("Prior mRNA","5.91x", "20.43x"),
  `BQ.1.1` = c("Prior AZ","0.78x", "7.54x"),
  ` ` = c("Prior mRNA","6.83x", "9.94x"),
  `XBB` = c("Prior AZ","0.30x", "3.27x"),
  ` ` = c("Prior mRNA","2.69x", "3.51x")
)

colnames(figure2_table) <- c("", "Ancestral/B.1"," ", "BA.1","  ", "BA.5","   ", "BQ.1.1","     " ,"XBB", "      ")

rownames(figure2_table) <- figure2_table[ ,1]

figure2_table<- subset(figure2_table, select=c("Ancestral/B.1"," ", "BA.1","  ", "BA.5","   ", "BQ.1.1","     " ,"XBB", "      "))

figure2 <- 
  ggplot(all_ic50_results4_b, aes(x=prime, y=logic50, color=serum_source))+
  facet_nested(.~pseudovirus,
               labeller = labeller(pseudovirus = pv.labs))+
  geom_boxplot(outlier.shape=NA, coef=Inf) + 
  theme_bw()+
  geom_point(position=position_jitterdodge(), alpha=0.3, size=1)+
  scale_color_manual(values=c("#000000","#ff9f49","#4c88ff"), labels=c("Healthy Controls", "HSCT (normal-responder)", "HSCT (low/non-responder)"))+
  theme(legend.position = "right",
        legend.title=element_blank(),
        plot.title=element_text(face="bold"))+
  scale_y_continuous(limits = c(0,6), breaks=c(0,1,2,3,4,5))+
  xlab(NULL)+
  ylab("log10 NT50")+
  geom_hline(yintercept=1.60206, linetype="dashed", color="black")+
  stat_pvalue_manual(stat.test_b, label="p.adj", size=3)

t1 <- ttheme_minimal(base_size = 9,
                     padding = unit(c(6, 3), "mm"),
                     core=list(
                       fg_params=list(fontface=c(rep("plain", 4))),
                       rowhead=list(),
                       bg_params = list(fill=c(NA, "#ff9f49","#4c88ff"), alpha=0.5)))

table_cowplot <- tableGrob(figure2_table[1:3, 1:10], theme=t1, rows=NULL)

ggdraw()+
  draw_plot(figure2, x=0,y=.2,width=1,height=.8)+
  draw_plot(table_cowplot, x=-0.33,y=0,width=1.5,height=.2)
#save as 11x8inch landscape

