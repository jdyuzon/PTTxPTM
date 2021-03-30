library(ggplot2)
library(ggpubr)
library(tidyr)

Input =("
SD      BGCs   gene
sig       11     20
nonsig     1      6
")

Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE,
                              row.names=1))
Matriz
fisher.test(Matriz,
            alternative="two.sided")
chisq.test(Matriz,
           correct=TRUE)
p.adjust(c(0.03866,0.3284,0.1974,0.6504,1.561e-05,0.000117,0.5224),
         method = "bonferroni")

##############################################################
##### Two Speed Genome #######################################
##############################################################
genomicfeatures<-read.table("/Users/jennifer/Desktop/Postdoc_Pteres/QuantGen/SegregationDistortion/genomicfeature.txt")
colnames(genomicfeatures)<-c("chrom","start","stop","genomicfeature")
genomicfeatures$length<-NULL
genomicfeatures$length[genomicfeatures$chrom=="chr1"]<-5917797
genomicfeatures$length[genomicfeatures$chrom=="chr2"]<-5608903
genomicfeatures$length[genomicfeatures$chrom=="chr3"]<-5204673
genomicfeatures$length[genomicfeatures$chrom=="chr4"]<-4483619
genomicfeatures$length[genomicfeatures$chrom=="chr5"]<-4379536
genomicfeatures$length[genomicfeatures$chrom=="chr6"]<-3318494
genomicfeatures$length[genomicfeatures$chrom=="chr7"]<-3273176
genomicfeatures$length[genomicfeatures$chrom=="chr8"]<-2553250
genomicfeatures$length[genomicfeatures$chrom=="chr9"]<-2371439
genomicfeatures$length[genomicfeatures$chrom=="chr10"]<-2252530
genomicfeatures$length[genomicfeatures$chrom=="chr11"]<-2056362
genomicfeatures$length[genomicfeatures$chrom=="chr12"]<-1287609
genomicfeatures$startdistance<-(genomicfeatures$start-0)/genomicfeatures$length
genomicfeatures$stopdistance<-(genomicfeatures$length-genomicfeatures$stop)/genomicfeatures$length
###minimum distance to telomere end scaled by length of chromosome
genomicfeatures$distance<-pmin(genomicfeatures$stopdistance,genomicfeatures$startdistance)*100
genomicfeatures$genomicfeature<-ordered(genomicfeatures$genomicfeature, levels = c("gene", "repeat", "BGC","effector"))

ggboxplot(subset(genomicfeatures, genomicfeature!="repeat"), x="genomicfeature",y="distance",
          fill="genomicfeature",palette =c("gray","deepskyblue3","palegreen3"))+
  stat_compare_means(comparisons = list(c("gene","BGC"),c("gene","effector")), 
                     method.args = list(alternative = "greater",method="wilcoxon.test"),size=2, label.x=1)+
  xlab("Genomic Feature")+ylab("Distance from Chromosome End (scaled by %percent length of chromosome)")+
  theme(legend.position = "none",
        axis.title=element_text(size=14,face="bold"))

p.val <- compare_means(distance ~ genomicfeature, alternative = "less",comparisons = list(c("gene","repeat"),c("gene","BGC"),c("gene","effector")), p.adjust.method = "bonferroni", method='wilcox.test', data = subset(genomicfeatures, genomicfeature!="repeat"))
p.adjust(c(2.22e-16,0.012,0.0045),
         method = "bonferroni")     ###6.66e-16 3.60e-02 1.35e-02

p.val <- p.val[1:2,] %>% mutate(y.position = c(53, 55.5))
Pathogenicity_Telomere<-ggboxplot(subset(genomicfeatures, genomicfeature!="repeat"), x="genomicfeature",y="distance",
          fill="genomicfeature",palette =c("gray","deepskyblue3","palegreen3"))+
  stat_pvalue_manual(p.val, label = "p.adj",size=2, label.x=1)+
  xlab("Genomic Feature")+ylab("Distance from Chromosome End
(scaled by %percent length of chromosome)")+
  theme(legend.position = "none",
        axis.title=element_text())


#### Segregation Distortion and the Two Speed Genome ########
df<-read.table("/Users/jennifer/Desktop/Postdoc_Pteres/QuantGen/SegregationDistortion/PTT_allgenomicregion_SD_density.txt")
colnames(df)<-c("chrom","start","stop","GeneDensity","NoPathGeneDensity","RepeatDensity","BGCsdensity","effectordensity", "SDdensity")
df$SD[df$SDdensity>0]<-"significant"
df$SD[df$SDdensity==0]<-"non-significant"
df$BGCs[df$BGCsdensity>0]<-"BGCs Region"
df$BGCs[df$BGCsdensity==0]<-"Genic Region"
df$effector[df$effectordensity>0]<-"Effector Region"
df$effector[df$effectordensity==0]<-"Genic Region"

###get distance of position from the telomere end
df$length<-NULL
df$length[df$chrom=="chr1"]<-5917797
df$length[df$chrom=="chr2"]<-5608903
df$length[df$chrom=="chr3"]<-5204673
df$length[df$chrom=="chr4"]<-4483619
df$length[df$chrom=="chr5"]<-4379536
df$length[df$chrom=="chr6"]<-3318494
df$length[df$chrom=="chr7"]<-3273176
df$length[df$chrom=="chr8"]<-2553250
df$length[df$chrom=="chr9"]<-2371439
df$length[df$chrom=="chr10"]<-2252530
df$length[df$chrom=="chr11"]<-2056362
df$length[df$chrom=="chr12"]<-1287609
df$startdistance<-(df$start-0)/df$length
df$stopdistance<-(df$length-df$stop)/df$length
###minimum distance to telomere end scaled by length of chromosome
df$distance<-pmin(df$stopdistance,df$startdistance)*100

df_long<-gather(data=df, key=genomicregion, value=density, c(GeneDensity:SDdensity), factor_key=TRUE)

## Segregation Distortion occurs at high repeat/low gene density
df_genrep<-subset(df_long, genomicregion=="NoPathGeneDensity"|genomicregion=="RepeatDensity"|genomicregion=="distance")

df_genrep$genomicregion<- factor(df_genrep$genomicregion, 
                               levels = c("GeneDensity","NoPathGeneDensity", "RepeatDensity","BGCsdensity","effectordensity","SDdensity","distance"),
                  labels = c("All Gene Density","Gene Density (per 50kb)","Repeat Density (per 50kb)","BGCsdensity","effectordensity","SDdensity","Distance from Chromosome End
(scaled by % length of chromosome)"))

SD_twospeed<-ggboxplot(df_genrep, x="SD",y="density",fill="SD",palette =c("blue","gray"))+
  facet_wrap(~genomicregion, scales = "free_y", strip.position = "left" )+
  stat_compare_means(size=2, label.x=1, label.y=100)+
  xlab("Segregation Distortion of 
       0-1xFGOB10Ptm-1 Cross (PTTxPTM)")+ ylab("")+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(),
        axis.title=element_text(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.key.size = unit(0.5,"cm"),
        legend.title = element_blank())


########### CDS: dnds pnps compared to segregation distortion
#MK_SD<-read.table("/Users/jennifer/Desktop/Postdoc_Pteres/Populations/MKtest/MK_results.pttptm.segregdistort.txt")
MK_SD<-read.table("/Users/jennifer/Desktop/Postdoc_Pteres/Populations/MKtest/MK_results.ptmptt.segregdistort.txt")
colnames(MK_SD)<-c("chrom","start","stop","CDS",
                   "ds","dn","ps","pn","dnds","pnps",
                   "NI","a","DoS",
                   "SD","allele")
MK_SD$allele[is.na(MK_SD$allele)]<-"non-significant"
MK_SD$allele<-factor(MK_SD$allele, levels = c("PTT_sd", "PTM_sd", "non-significant"))
MK_SD$mk<-MK_SD$dnds/MK_SD$pnps

mean(na.omit(MK_SD$dnds)) #PTT: 0.1973684 #PTM: 0.1973684
mean(na.omit(MK_SD$pnps)) #PTT: 0.8442814 #PTM: 0.8442814
mean(na.omit(MK_SD$mk[MK_SD$mk!="Inf"])) #PTT: 0.5104612 #PTM: 0.5104612

median(na.omit(MK_SD$dnds)) #PTT: 0 #PTM: 0
median(na.omit(MK_SD$pnps)) #PTT: 0.5333333 #PTM: 0.5333333
median(na.omit(MK_SD$mk[MK_SD$mk!="Inf"])) #PTT: 0 #PTM: 0

sum(na.omit(MK_SD$dn)) #PTT: 115  #PTM: 115
sum(na.omit(MK_SD$ds)) #PTT: 151  #PTM: 151
sum(na.omit(MK_SD$pn)) #PTT: 20327 #PTM: 20327
sum(na.omit(MK_SD$ps)) #PTT: 33408 #PTM: 33408

hist<-list()
for (stat in c("dnds","pnps","NI","a","DoS")){
hist[[stat]]<-ggplot(MK_SD, aes_string(x=stat, fill="allele")) +
  geom_histogram(position="dodge")+
  theme(legend.position="top")+
  scale_fill_manual(values=c("red","blue","gray"))+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(),
        axis.title=element_text(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.key.size = unit(0.5,"cm"),
        legend.title = element_blank())
}

boxplot<-list()
for (stat in c("dnds","pnps","NI","a","DoS","mk")){
boxplot[[stat]]<-ggboxplot(MK_SD, x="SD",y=stat,fill="SD",palette =c("blue","gray"))+
  stat_compare_means(size=2,
                     method="wilcox.test", method.args = list(alternative = "greater"),
                     comparisons = list(c("significant","non-significant")))+
  xlab("")+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        #strip.text = element_text(),
        axis.title=element_text(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.key.size = unit(0.5,"cm"),
        legend.title = element_blank())
}

p<-ggarrange(boxplot[["dnds"]],boxplot[["pnps"]],
          legend="bottom", common.legend = TRUE)

pdf("~/Desktop/Postdoc_Pteres/Populations/MKtest/Figure_SegregDist_dnds_pnps.pdf",
    width=5, height=5)
plot(p)
dev.off()

## Towards the beginning and ends of chromsomes high repeat/low gene density
Telomere_GeneDensity<-ggplot(df, aes(x=distance, y=NoPathGeneDensity,)) +  
  geom_point(color="gray", alpha=0.4)+geom_smooth(se = FALSE,colour = "red")+
  #xlab("Distance from Telomere (scaled by % length of chromosome)")+
  annotate("text", x=10, y=100, label= "tau = 0.15  p.adj = 8.80e-16", size=3)+
  xlab("")+ ylab("Gene Density (per 50kb)")+
  theme_classic()
Telomere_RepeatDensity<-ggplot(df, aes(x=distance, y=RepeatDensity)) +  
  geom_point(color="gray", alpha=0.4)+geom_smooth(se = FALSE,colour = "red")+
  annotate("text", x=10, y=70, label= "tau = -0.18  p.adj = 8.80e-16", size=3)+
  xlab("Distance from Telomere (scaled by % length of chromosome)")+
  ylab("Repeat Density (per 50kb)")+
  theme_classic()

cor.test( ~ NoPathGeneDensity + distance,data=df,method = "kendall") #p-value < 2.2e-16     p.adj= 8.8000e-16*   tau 0.1510859
cor.test( ~ RepeatDensity + distance,data=df,method = "kendall")     #p-value < 2.2e-16     p.adj= 8.8000e-16*   tau -0.1829826 
cor.test( ~ BGCsdensity + distance,data=df,method = "kendall")       #p-value = 1.818e-06   p.adj= 7.2720e-06*   tau -0.06130693
cor.test( ~ effectordensity + distance,data=df,method = "kendall")   #p-value = 0.03713     p.adj= 1.4852e-01    tau -0.02653613
p.adjust(c(2.2e-16, 2.2e-16, 1.818e-06, 0.03713), method = "bonferroni")

SD_twospeed_plots<-ggarrange(Telomere_GeneDensity,
          Pathogenicity_Telomere,
          Telomere_RepeatDensity,
          SD_twospeed,
          widths = c(1.7,2))

pdf("~/Desktop/Postdoc_Pteres/QuantGen/SegregationDistortion/Figure_SegregDist_twospeed.pdf",
    width=11.9, height=14)
SD_twospeed_plots
dev.off()



## repeat and gene density are inversely related
ggplot(df, aes(x=NoPathGeneDensity, y=RepeatDensity, color=start)) +
  geom_point()+geom_smooth(method=lm)+facet_wrap(~chrom)


