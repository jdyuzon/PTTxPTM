rm(list = ls())
library(tidyverse)
library(ggpubr)
require(latex2exp)

##### Genomic Features
chrom_names<-c("1","2","3","4","5","6","7","8","9","10","11","12","m")
goodChrOrder <- c(paste(c(1:12),sep=""),"mitochondria")
names(chrom_names)<-goodChrOrder

genedensity<-read.table("/home/yuzon/PTTxPTMtest/PTTxPTMcompare/SegregationDistortion_files/PTTxPTM_tracks.genedensity.bed", header = FALSE,sep = "\t")
colnames(genedensity)<-c("chrom","start","stop","genecount")
genedensity<-genedensity[grepl("chrPTT|m86|contig", genedensity$chrom), , drop = FALSE]
genedensity<-data.frame(lapply(genedensity, function(x) {gsub("chrPTT_", "", x)}))
genedensity<- data.frame(lapply(genedensity, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
genedensity[,2:4] <- data.frame(lapply(genedensity[,2:4], function(x) {as.numeric(as.character(x))}))
genedensity_genome<-genedensity[grepl(c("1|2|3|4|5|6|7|8|9|10|11|12|mitochondria$"),genedensity$chrom), , drop = FALSE]
genedensity_contig<-genedensity[grepl("0-1_contig", genedensity$chrom), , drop = FALSE]
genedensity_genome$chrom <- factor(genedensity_genome$chrom ,levels=goodChrOrder)
genedensity_genome$stat<-"Fst"  #####add track to last facet of the divergence graph
genedensity_genome$stat <- factor(genedensity_genome$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst","ABBA","BABA","D","fd","fdM"))

####format effectors
effectors<-read.table("/home/yuzon/PTTxPTMtest/PTTxPTMcompare/SegregationDistortion_files/PTTxPTM_tracks.effectors.bed", header = FALSE,sep = "\t")
colnames(effectors)<-c("chrom","start","stop")
effectors<-effectors[grepl("chrPTT|m86|contig", effectors$chrom), , drop = FALSE]
effectors<-data.frame(lapply(effectors, function(x) {gsub("chrPTT_", "", x)}))
effectors<- data.frame(lapply(effectors, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
effectors[,2:3] <- data.frame(lapply(effectors[,2:3], function(x) {as.numeric(as.character(x))}))
effectors$type<-"effector"
effectors_genome<-effectors[grepl(c("1|2|3|4|5|6|7|8|9|10|11|12|mitochondria$"), effectors$chrom), , drop = FALSE]
effectors_contig<-effectors[grepl("0-1_contig", effectors$chrom), , drop = FALSE]
effectors_genome$chrom <- factor(effectors_genome$chrom ,levels=goodChrOrder)
effectors_genome$stat<-"Fst"
effectors_genome$stat <- factor(effectors_genome$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst","ABBA","BABA","D","fd","fdM"))

####format biosynthetic clusters
bsc<-read.table("/home/yuzon/PTTxPTMtest/PTTxPTMcompare/SegregationDistortion_files/PTTxPTM_tracks.bsc.bed", header = FALSE,sep = "\t")
colnames(bsc)<-c("chrom","start","stop","region")
bsc<-bsc[grepl("chrPTT|m86|contig", bsc$chrom), , drop = FALSE]
bsc<-data.frame(lapply(bsc, function(x) {gsub("chrPTT_", "", x)}))
bsc<- data.frame(lapply(bsc, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
bsc[,2:3] <- data.frame(lapply(bsc[,2:3], function(x) {as.numeric(as.character(x))}))
bsc$types<-"bsc"
bsc_genome<-bsc[grepl(c("1|2|3|4|5|6|7|8|9|10|11|12|mitochondria$"), bsc$chrom), , drop = FALSE]
bsc_contig<-bsc[grepl("0-1_contig", bsc$chrom), , drop = FALSE]
bsc_genome$chrom <- factor(bsc_genome$chrom ,levels=goodChrOrder)
bsc_genome$stat<-"Fst"
bsc_genome$stat <- factor(bsc_genome$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst","ABBA","BABA","D","fd","fdM"))

####format mat loci
mat<-read.table("/home/yuzon/PTTxPTMtest/PTTxPTMcompare/SegregationDistortion_files/matloci.bed", header = FALSE,sep = "\t")
colnames(mat)<-c("chrom","start","stop")
mat<-mat[grepl("chrPTT|m86|contig", mat$chrom), , drop = FALSE]
mat<-data.frame(lapply(mat, function(x) {gsub("chrPTT_", "", x)}))
mat[,2:3] <- data.frame(lapply(mat[,2:3], function(x) {as.numeric(as.character(x))}))
mat$types<-"mat"
mat$chrom <- factor(mat$chrom ,levels=goodChrOrder)
mat$stat<-"Fst"
mat$stat <- factor(mat$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst","ABBA","BABA","D","fd","fdM"))

########################################################
####https://github.com/simonhmartin/genomics_general####
####http://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/divergence-and-gene-flow-exercise/
########################################################
div_data <- read.csv("divergence.csv.tmp", as.is=T)
div_data$scaffold<-gsub("chr","",div_data$scaffold)
div_data$scaffold <- factor(div_data$scaffold,levels=goodChrOrder)
div_df<-reshape(div_data, 
                     direction = "long",
                     varying = list(names(div_data)[6:9]),
                     v.names = "value",
                     timevar = "stat",
                     times = c("pi_ptt","pi_ptm","dxy_ptt_ptm","Fst_ptt_ptm"),
                     idvar = c("scaffold", "start","end","mid","sites"))


g<-ggplot() +
  geom_line(data=div_df, aes(mid/10^6, value, colour = stat))+ 
  theme_light() +   
  facet_grid(stat~scaffold,scales="free",labeller = labeller(.rows = label_parsed))+
  theme(legend.position = "none")

### Get the Dxy and Fst calculation for the whole genome
div_df$stat <- factor(div_df$stat, levels = c("pi_ptm", "pi_ptt", "dxy_ptt_ptm","Fst_ptt_ptm"), 
                      labels = c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
div_df$stat <- factor(div_df$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
div_df_ptm <- subset(div_df,stat==as.character(expression(paste("PTM ",pi))))
div_df_ptt <- subset(div_df,stat==as.character(expression(paste("PTT ",pi))))
div_df_dxy <- subset(div_df,stat=='Dxy')
div_df_fst <- subset(div_df,stat=='Fst')
div_df_ptm$stat <- factor(div_df_ptm$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
div_df_ptt$stat <- factor(div_df_ptt$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
div_df_dxy$stat <- factor(div_df_dxy$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
div_df_fst$stat <- factor(div_df_fst$stat, c(expression(paste("PTM ", pi)), expression(paste("PTT ", pi)), "Dxy","Fst"))
### Identify Outliers as the upper 0.1% quantiles and get these regions
options(scipen=999)

outlier_ptm_cutoff<-quantile(na.omit(div_df_ptm$value), 0.90)
outlier_ptm_regions<-subset(div_df_ptm,value>=outlier_ptm_cutoff)

outlier_ptt_cutoff<-quantile(na.omit(div_df_ptt$value), 0.90)
outlier_ptt_regions<-subset(div_df_ptt,value>=outlier_ptt_cutoff)

outlier_dxy_cutoff<-quantile(na.omit(div_df_dxy$value), 0.90)
outlier_dxy_regions<-subset(div_df_dxy,value>=outlier_dxy_cutoff)
write.table(outlier_dxy_regions[,c("scaffold","start","end")], file = "outlier_dxy_regions.bed", sep = "\t",row.names = FALSE, col.names = FALSE,quote=FALSE)

outlier_fst_cutoff<-quantile(na.omit(div_df_fst$value), 0.90)
outlier_fst_regions<-subset(div_df_fst,value>=outlier_fst_cutoff)
write.table(outlier_fst_regions[,c("scaffold","start","end")], file = "outlier_fst_regions.bed", sep = "\t",row.names = FALSE, col.names = FALSE,quote=FALSE)

pdf("DivergenceStats.pdf", width=24, height=8)
plot(g)
graphics.off()

################################################################
colnames(genedensity_genome)<-c("scaffold","start", "end","genecount","stat")
colnames(effectors_genome)<-c("scaffold","start", "end","types","stat")
colnames(bsc_genome)<-c("scaffold", "start","end","region","types","stat")
colnames(mat)<-c("scaffold","start","end","types","stat")

divergence_plot<-ggplot() +
  geom_line(data=div_df, aes(mid/10^6, value, colour = stat))+ 
  theme_light() +   
  facet_grid(stat~scaffold,scales="free",labeller = labeller(.rows = label_parsed))+
  theme(legend.position = "none")+
  geom_hline(data=div_df_ptm, aes(yintercept=outlier_ptm_cutoff), linetype="dashed", color = "black")+
  geom_hline(data=div_df_ptt, aes(yintercept=outlier_ptt_cutoff), linetype="dashed", color = "black")+
  geom_hline(data=div_df_fst, aes(yintercept=outlier_fst_cutoff), linetype="dashed", color = "black")+
  geom_hline(data=div_df_dxy, aes(yintercept=outlier_dxy_cutoff), linetype="dashed", color = "black")+
  geom_hline(data=div_df_fst, aes(yintercept = -0.35),color="white",size = 25)+
  xlab("Nuclear Genome (position Mb)")+
  theme_light()+ theme(panel.margin = unit(0, "lines"),
                      panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
                      #strip.text.x.top = element_text(angle = 90),
                      legend.position = "none",
                      strip.text.x.top = element_text(size = 15, color = "black",face = "bold"),
                      axis.ticks.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title=element_text(size=12,face="bold"))+
  geom_rect(data=genedensity_genome[ which(genedensity_genome$scaffold !='mitochondria' ),],aes(xmin = start/10^6, xmax = end/10^6,ymin=-0.05,
                                                                                             ymax=-0.1,fill=genecount),show.legend = FALSE)+ 
  scale_fill_gradient(low="#FFFFFF00",high="gray19")+
  geom_rect(data=effectors_genome[ which(effectors_genome$scaffold !='mitochondria' ),],aes(xmin = start/10^6, xmax = end/10^6,ymin=-0.15,
                                                                                         ymax=-0.2, colour=types),show.legend = FALSE)+
  geom_linerange(data=bsc_genome[ which(bsc_genome$scaffold !='mitochondria' ),],aes(xmin = start/10^6, xmax = end/10^6,y=-0.3,
                                                                                  colour=types),size=3.5,show.legend = FALSE)+
  geom_linerange(data=mat[ which(mat$scaffold !='mitochondria' ),],aes(xmin = (start-10000)/10^6, xmax = (end+10000)/10^6,y=-0.4,
                                                                    colour=types),size=5,alpha=1,show.legend = FALSE)+
  scale_colour_manual(limits = c(expression(paste("PTT ", pi)), expression(paste("PTM ", pi)), "Dxy","Fst","bsc","effector","mat"), 
                      values=c("red", "blue", "darkgreen","orange","deepskyblue3","palegreen3","purple"))
divergence_plot

g <- ggplot_gtable(ggplot_build(divergence_plot))
strips <- which(grepl('strip-', g$layout$name))
pal <- rev(c( "violet","blue","darkgreen","orange", "red", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2"))
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  #l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}


plot(g)

pdf("DivergenceStats.shm.pdf", width=24, height=8)
plot(g)
graphics.off()
