library(ggpubr)
library(tidyverse)
library(grid)
library(car)

phenotypes<-read.csv("PterespopsonKombar.csv", header = TRUE)
colnames(phenotypes)[8]<-"mean.percent.leaf.coverage"
phenotypes$labels <- paste0(phenotypes$CrossParents)
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='N']<-'Net Lesion'
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='S']<-'Spot Lesion'
phenotypes$Net.blotch.form <- factor(phenotypes$Net.blotch.form,
                                levels = c('Net Lesion','Spot Lesion'),ordered = TRUE)

phenotypes$Cross.Type <- factor(phenotypes$Cross.Type,
                                levels = c('Parents','Intra','Inter'),ordered = TRUE)

phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"]<-phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"]/max(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"])
phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"]<-phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"]/max(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"])



phenotypes$mean.percent.leaf.coverage<-phenotypes$mean.percent.leaf.coverage/max(na.omit(phenotypes$mean.percent.leaf.coverage))
phenotypes$RT.Average<-phenotypes$RT.Average/max(na.omit(phenotypes$RT.Average))

my_comparisons <- list(c("Intra", "Inter"), c("Parents", "Inter"), c("Intra", "Parents"))

a<-ggboxplot(phenotypes, x = "Cross.Type", y = "percent.leaf.coverage", add = "jitter")+ 
  facet_wrap(~Rep)+
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  xlab("")+ylab("Fitness (% Leaf Coverage)")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) 

A<-ggboxplot(phenotypes, x = "Cross.Type", y = "Kombar.reaction.type", add = "jitter")+ 
  facet_wrap(~Rep)+
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  xlab("")+ylab("Fitness (Disease Score)")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))



b<-ggboxplot(phenotypes, x = "Cross.Type", y = "percent.leaf.coverage", 
             color = "Net.blotch.form",palette =c("red", "blue"),add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons , size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  facet_grid(Rep~Net.blotch.form, scales = "free_x")+
  xlab("")+ylab("Fitness (% Leaf Coverage)")+
  theme(legend.position="None", strip.text.x = element_text(color = "white"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))


p <- ggplot_gtable(ggplot_build(b))
strips <- which(grepl('strip-', p$layout$name))

pal<-c("red","blue","grey53", "grey74","grey94")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', p$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  p$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

plot(p)



pdf("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/Phenotype_percentleaf_allrep")
plot(p)
dev.off()


####Medians-percent.leaf.coverage
###Parents-Net
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
             phenotypes$Rep==1 &
             phenotypes$Cross.Type=='Parents']) #0.875
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Parents']) #0.875
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Parents']) #0.6

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                       phenotypes$Cross.Type=='Parents'])) #0.7263158

###Intraspecies-Net
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Intra']) #0.7
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Intra']) #0.6
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Intra']) #0.7

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                                       phenotypes$Cross.Type=='Intra'])) #0.6315789

###Interspecies-Net
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Inter']) #0.2
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Inter']) #0.35
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Inter']) #0.25

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                                       phenotypes$Cross.Type=='Inter'])) #0.3473684

###Parents-Spot
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Parents']) #0.325
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Parents']) #0.4
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Parents']) #0.65

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Parents'])) #0.4631579

###Intraspecies-Spot
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Intra']) #0.4
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Intra']) #0.4
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Intra']) #0.5

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Intra'])) #0.4421053

###Interspecies-Spot
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Inter']) #0.2
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Inter']) #0.225
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Inter']) #0.2

median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Inter'])) #0.2736842


leveneTest(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#Df F value Pr(>F)
#group   2   0.766 0.4657
#334   
fligner.test(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#data:  percent.leaf.coverage by as.factor(Rep)
#Fligner-Killeen:med chi-squared = 1.2318, df = 2, p-value = 0.5402

###################################################################
##################### Control for Lesion Type #####################
###################################################################
phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"]<-phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"]/max(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Spot Lesion"]))
phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"]<-phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"]/max(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=="Net Lesion"]))



####Pooled data (all reps)
b<-ggboxplot(phenotypes, x = "Cross.Type", y = "mean.percent.leaf.coverage", 
             color = "Net.blotch.form",palette =c("red", "blue"),add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons , size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  facet_grid(~Net.blotch.form, scales = "free_x")+
  xlab("")+ylab("Fitness (% Leaf Coverage)")+
  theme(legend.position="None", strip.text.x = element_text(color = "white"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))


p <- ggplot_gtable(ggplot_build(b))
strips <- which(grepl('strip-', p$layout$name))

pal<-c("red","blue","grey53", "grey74","grey94")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', p$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  p$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

plot(p)


pdf("Phenotype_percentleaf_pooled")
plot(p)
dev.off()


