library(ggpubr)
library(tidyverse)
library(grid)
library(car)

phenotypes<-read.csv("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/PterespopsonKombar.csv", header = TRUE)
colnames(phenotypes)[8]<-"mean.percent.leaf.coverage"
phenotypes$labels <- paste0(phenotypes$CrossParents)
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='N']<-'Net Lesion'
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='S']<-'Spot Lesion'
phenotypes$Net.blotch.form <- factor(phenotypes$Net.blotch.form,
                                levels = c('Net Lesion','Spot Lesion'),ordered = TRUE)

phenotypes$Cross.Type <- factor(phenotypes$Cross.Type,
                                levels = c('Parents','Intra','Inter'),ordered = TRUE)

phenotypes$percent.leaf.coverage<-phenotypes$percent.leaf.coverage/max(phenotypes$percent.leaf.coverage)
phenotypes$Kombar.reaction.type<-phenotypes$Kombar.reaction.type/max(phenotypes$Kombar.reaction.type)

phenotypes$mean.percent.leaf.coverage<-phenotypes$mean.percent.leaf.coverage/max(na.omit(phenotypes$mean.percent.leaf.coverage))
phenotypes$RT.Average<-phenotypes$RT.Average/max(na.omit(phenotypes$RT.Average))

my_comparisons <- list(c("Intra", "Inter"), c("Parents", "Inter"), c("Intra", "Parents"))

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



B<-ggboxplot(phenotypes, x = "Cross.Type", y = "Kombar.reaction.type", 
             color = "Net.blotch.form",palette =c("red", "blue"),add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  facet_grid(Rep~Net.blotch.form, scales = "free_x")+
  xlab("")+ylab("Fitness (Disease Score)")+
  theme(legend.position="None", strip.text.x = element_text(color = "white"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))  


g <- ggplot_gtable(ggplot_build(B))
strips <- which(grepl('strip-', g$layout$name))

pal<-c("red","blue","grey53", "grey74","grey94")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

plot(g)


pdf("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/Phenotype_percentleaf_allrep")
plot(p)
dev.off()

pdf("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/Phenotype_DiseaseScore_allrep")
plot(g)
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


####Medians-DiseaseScore
###Parents-Net
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Parents']) #0.9
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Parents']) #0.825
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Parents']) #0.75

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Net Lesion' &
                                         phenotypes$Cross.Type=='Parents'])) #0.8660236

###Intraspecies-Net
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Intra']) #0.75
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Intra']) #0.7
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Intra']) #0.75

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Net Lesion' &
                                       phenotypes$Cross.Type=='Intra'])) #0.7856377


###Interspecies-Net
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Inter']) #0.35
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Inter']) #0.45
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Net Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Inter']) #0.4

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Net Lesion' &
                                       phenotypes$Cross.Type=='Inter'])) #0.4287245

###Parents-Spot
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Parents']) #0.35
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Parents']) #0.35
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Parents']) #0.45

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Spot Lesion' &
                                       phenotypes$Cross.Type=='Parents'])) #0.375134

###Intraspecies-Spot
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Intra']) #0.35
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Intra']) #0.35
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Intra']) #0.35

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Spot Lesion' &
                                       phenotypes$Cross.Type=='Intra'])) #0.375134

###Interspecies-Spot
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==1 &
                                          phenotypes$Cross.Type=='Inter']) #0.25
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==2 &
                                          phenotypes$Cross.Type=='Inter']) #0.225
median(phenotypes$Kombar.reaction.type[phenotypes$Net.blotch.form=='Spot Lesion' &
                                          phenotypes$Rep==3 &
                                          phenotypes$Cross.Type=='Inter']) #0.225

median(na.omit(phenotypes$RT.Average[phenotypes$Net.blotch.form=='Spot Lesion' &
                                       phenotypes$Cross.Type=='Inter'])) #0.249732


#### Test for homogeneity
leveneTest(Kombar.reaction.type ~ as.factor(Rep), data = phenotypes)
#Df F value Pr(>F)
#group   2  0.9944  0.371
#334 
fligner.test(Kombar.reaction.type ~ as.factor(Rep), data = phenotypes)
#data:  Kombar.reaction.type by as.factor(Rep)
#Fligner-Killeen:med chi-squared = 1.9351, df = 2, p-value = 0.38

leveneTest(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#Df F value Pr(>F)
#group   2   0.766 0.4657
#334   
fligner.test(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#data:  percent.leaf.coverage by as.factor(Rep)
#Fligner-Killeen:med chi-squared = 1.2318, df = 2, p-value = 0.5402

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



B<-ggboxplot(phenotypes, x = "Cross.Type", y = "RT.Average", 
             color = "Net.blotch.form",palette =c("red", "blue"),add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  facet_grid(~Net.blotch.form, scales = "free_x")+
  xlab("")+ylab("Fitness (Disease Score)")+
  theme(legend.position="None", strip.text.x = element_text(color = "white"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))  


g <- ggplot_gtable(ggplot_build(B))
strips <- which(grepl('strip-', g$layout$name))

pal<-c("red","blue","grey53", "grey74","grey94")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

plot(g)

pdf("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/Phenotype_percentleaf_pooled")
plot(p)
dev.off()

pdf("~/Desktop/Postdoc_Pteres/QuantGen/Pteres_Postzygotic_F1_phenotypes/FinalPhenotypingResults/Phenotype_DiseaseScore_pooled")
plot(g)
dev.off()

