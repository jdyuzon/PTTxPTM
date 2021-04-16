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


pdf("Phenotype_percentleaf_allrep")
plot(p)
dev.off()

pdf("Phenotype_DiseaseScore_allrep")
plot(g)
dev.off()


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

pdf("Phenotype_percentleaf_pooled")
plot(p)
dev.off()

pdf("Phenotype_DiseaseScore_pooled")
plot(g)
dev.off()

