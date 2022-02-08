library(vcfR)
library(onemap)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyr)
library(LDheatmap)
library(RColorBrewer)

############################################################################
#####################Recombination Blocks: Individuals######################
############################################################################

############ Mapping to PTT and PTM genomes ##############################
goodChrOrder <- c(paste("chr",c(1:12),sep=""),"mitochondria")
goodIsoOrder <- c(1:70,72:79,"79b",80:120)


#### Finer Quality (extend ranges to nearest similar allele): PTM coordinates mapped to PTT 
blockbed<-read.table("all.common_unique.blocks.bed.tmp", header = FALSE,sep = "\t",fill=TRUE)

#### isolate order for publication
goodIsoOrder <- c("1","2","3","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","25","27","28","29","30","31","32","33","34","36","37","38","39","40","41","42","43","45","46","47","48","49","50","51","52","53","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","72","75","76","77","79b","81","82","83","84","85","86","87","88","90","91","92","93","95","96","97","98","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120")

blockbed$sizes<-blockbed$V3-blockbed$V2
blockbed$V1 <- factor(blockbed$V1,levels=goodChrOrder)
blockbed$V5 <- factor(blockbed$V5,levels=goodIsoOrder)

ggplot()+
  geom_rect(data=na.omit(blockbed), mapping=aes(xmin=V2, xmax=V3, ymin=0, ymax=1, fill=V4), alpha=0.5)+
  scale_fill_manual(breaks = c("PTT", "PTM","common"),values=c("red", "blue","purple"))+
  facet_grid(V5~V1,scales="free", space="free_x", switch="both")+
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=8, angle=90),
        strip.text.y.left = element_text(angle = 0),
        panel.spacing=unit(0, "lines"))


#####Histogram of block sizes by chromosomes
blockbed$sizes<-blockbed$V3-blockbed$V2
ggplot(blockbed, aes(x=log10(sizes), fill=V4)) +
  geom_histogram(position="identity", alpha=0.5)+
  scale_fill_manual(breaks = c("PTT", "PTM","common"),values=c("red", "blue","purple"))+
  facet_wrap(~V1,scales="free")

####################### Mapping to PTT genome ##############################

#### Finer Quality (extend ranges to nearest similar allele): PTM coordinates mapped to PTT 
blockPTT<-read.table("all.common_unique.blocks.bed", header = FALSE,sep = "\t")
blockPTT_matmito<-read.table("all.common_unique.blocks.matmito.bed", header = FALSE,sep = "\t")

blockPTT$V1 <- factor(blockPTT$V1,levels=goodChrOrder)
blockPTT$V5 <- factor(blockPTT$V5,levels=goodIsoOrder)
blockPTT_matmito$V1 <- factor(blockPTT_matmito$V1,levels=goodChrOrder)
blockPTT_matmito$V5 <- factor(blockPTT_matmito$V5,levels=goodIsoOrder)


ggplot()+
  geom_rect(data=na.omit(blockPTT), mapping=aes(xmin=V2, xmax=V3, ymin=0, ymax=1, fill=V4), alpha=0.5)+
  scale_fill_manual(breaks = c("PTT", "PTM","common"),values=c("red", "blue","purple"))+
  facet_grid(V5~V1,scales="free", space="free_x", switch="both")+
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=8, angle=90),
        strip.text.y.left = element_text(angle = 0),
        panel.spacing=unit(0, "lines"))


ggplot()+
  geom_rect(data=blockPTT_matmito, mapping=aes(xmin=V2, xmax=V3, ymin=0, ymax=1, fill=V4), alpha=0.5)+
  scale_fill_manual(breaks = c("PTT", "PTM","common"),values=c("red", "blue","purple"))+
  facet_grid(V5~V1,scales="free",  switch="both")+
  theme(axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=8, angle=0),
        strip.text.y.left = element_text(angle = 0),
        panel.spacing=unit(0, "lines"))

############################################################################
######################Segregation Distortion: Hybrid Population#############
############################################################################
#df<-read.table("~/Downloads/all.recomblocks.bed.tmp", header = TRUE,sep = "\t",na.strings = c(".","0"))
df<-read.table("all.recomblocks.bed.tmp", header = TRUE,sep = "\t",na.strings = c(".","0")) #,fill=TRUE
df$chrom <- factor(df$chrom,levels=c(goodChrOrder,
                                     "0-1_contig_32","0-1_contig_39",
                                     "0-1_contig_40","0-1_contig_43",
                                     "0-1_contig_44","0-1_contig_47",
                                     "0-1_contig_48","0-1_contig_50",
                                     "0-1_contig_51","0-1_contig_52",
                                     "0-1_contig_53","0-1_contig_54",
                                     "0-1_contig_55","0-1_contig_56",
                                     "0-1_contig_57","0-1_contig_59",
                                     "0-1_contig_60","0-1_contig_61",
                                     "0-1_contig_62","0-1_contig_63",
                                     "0-1_contig_65","0-1_contig_67",
                                     "0-1_contig_68","0-1_contig_71",
                                     "0-1_contig_72","0-1_contig_73",
                                     "0-1_contig_74","0-1_contig_75",
                                     "0-1_contig_78","0-1_contig_79",
                                     "0-1_contig_80","0-1_contig_81",
                                     "0-1_contig_82","0-1_contig_84"
))

df_PTT<-df[c(1,2,3,grep("PTT.support", names(df)))]
df_PTM<-df[c(1,2,3,grep("PTM.support", names(df)))]

df$obs_ptt<-apply(df_PTT[3:109], 1, function(x) length(which(x=="PTT")))
df$obs_ptm<-apply(df_PTM[3:109], 1, function(x) length(which(x=="PTM")))
df$sum<-df$obs_ptt+df$obs_ptm
df <- subset(df, sum>=90)
df$exp_ptt<-0.5*df$sum
df$exp_ptm<-0.5*df$sum
df$freq_ptt<-df$obs_ptt/df$sum
df$freq_ptm<-df$obs_ptm/df$sum
df$prob_ptt<-0.5
df$prob_ptm<-0.5
df$x2<-((df$obs_ptt-df$exp_ptt)/df$exp_ptt)^2+((df$obs_ptm-df$exp_ptm)/df$exp_ptm)^2

df1<-data.frame(df$obs_ptt, df$obs_ptm, df$prob_ptt, df$prob_ptm)

get_chisq <- function(x, prbs) {
  chsq <- chisq.test(x=x, p=prbs)
  ans <- cbind(statistic=chsq$statistic[[1]],
               df=chsq$parameter[[1]],
               p.value=chsq$p.value)
  ans
}

sol<-data.frame(t(apply(df1, 1, function(x) get_chisq(x[1:2], x[3:4]))))
names(sol)<-c("statistic","df","p.value")
sol$Bonferroni =
  p.adjust(sol$p.value,
           method = "bonferroni")
segreg_test<-cbind(sol$statistic,sol$p.value,df$sum/10*100)

segreg_test<-data.frame(apply(segreg_test, 2, function(x) as.numeric(as.character(x))))
segreg_test<-cbind(df$chrom,df$start,df$stop,rownames(df),df$freq_ptt, df$freq_ptm,"1:1", segreg_test)
colnames(segreg_test)<-c("chrom","start","stop","Marker","freq_ptt","freq_ptm","H0","Chi-square","p-value", "% genot.")
segreg_test$bonferonni<-sol$Bonferroni
segreg_test$significant<-ifelse(sol$Bonferroni<0.05,"significant","non-significant")

segreg_test_table<-segreg_test
segreg_test_table$SD_parent[segreg_test_table$freq_ptt>0.5 & segreg_test_table$significant=="significant"]<-"PTT_sd"
segreg_test_table$SD_parent[segreg_test_table$freq_ptt<0.5 & segreg_test_table$significant=="significant"]<-"PTM_sd"
segreg_test_table$SD_parent[segreg_test_table$significant=="non-significant"]<-"NA"
#write.table(segreg_test_table, file = "~/Desktop/Postdoc_Pteres/QuantGen/SegregationDistortion/all.segregdistort.PTTxPTM.txt", append = FALSE, quote = FALSE, sep = "\t",
#            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
#            col.names = FALSE, qmethod = c("escape", "double"),
#            fileEncoding = "")

segreg_genome<-segreg_test[grepl("chr|mito", segreg_test$chrom), , drop = FALSE]
segreg_contig<-segreg_test[grepl("contig", segreg_test$chrom), , drop = FALSE]

ptm_seggenome<-ggplot(segreg_genome, aes(y=freq_ptm)) +
  geom_linerange(aes(xmin = start, xmax = stop, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("gray", "blue"))+
  facet_grid(~chrom, scales = "free_x", space = "free_x")+
  geom_hline(yintercept = 0.5)+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))

ptt_seggenome<-ggplot(segreg_genome, aes(y=freq_ptt)) +
  geom_linerange(aes(xmin = start, xmax = stop, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("gray", "red"))+
  facet_grid(~chrom, scales = "free_x", space = "free_x")+
  geom_hline(yintercept = 0.5)+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))

ptm_segcontig<-ggplot(segreg_contig, aes(y=freq_ptm)) +
  geom_linerange(aes(xmin = start, xmax = stop, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("gray", "blue"))+
  facet_grid(~chrom, scales = "free_x", space = "free_x")+
  geom_hline(yintercept = 0.5)+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))

ptt_segcontig<-ggplot(segreg_contig, aes(y=freq_ptt)) +
  geom_linerange(aes(xmin = start, xmax = stop, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("gray", "red"))+
  facet_grid(~chrom, scales = "free_x", space = "free_x")+
  geom_hline(yintercept = 0.5)+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggarrange(ptt_seggenome,ptm_seggenome,ptt_segcontig,ptm_segcontig, 
          ncol = 1, nrow = 4)


############################################################################
###################### Recombination Hotspots ##############################
############################################################################

rbreak_dens<-read.table("all.breaks.density", header = FALSE,sep = "\t")
##rbreak<-read.table("~/Downloads/testing_ground/all.breaks.blocks", header = FALSE,sep = "\t")
###rbreak_dens<-read.table("~/Downloads/all.breaks.density", header = FALSE,sep = "\t")

rbreak<-rbreak_dens

rbreak_dens <- subset(rbreak_dens,V4>0)

rbreak_dens<-(unique(rbreak_dens[,1:4]))
median(rbreak_dens$V4)
mean(rbreak_dens$V4)
#https://www.genetics.org/content/201/3/1213#sec-1
#Poisson probability: lambda = mean occurance per interval
#http://www.r-tutor.com/elementary-statistics/probability-distributions/poisson-distribution
ppois(7, lambda=mean(rbreak_dens$V4),lower=FALSE)
rbreak_dens$ppois<-ppois(rbreak_dens$V4, lambda=mean(rbreak_dens$V4),lower=FALSE)
rbreak_dens$Bonferroni=p.adjust(rbreak_dens$ppois,method = "bonferroni")
rbreak_dens$significant<-ifelse(rbreak_dens$Bonferroni<=0.05, "hotspot", "non-hotspot")

rbreak<-merge(rbreak_dens,rbreak)
rbreak<-unique(rbreak[,7:13])

rbreak$size<-rbreak$V7-rbreak$V6
mu <- ddply(rbreak, "significant", summarise, grp.mean=mean(V10/1000))
med <- ddply(rbreak, "significant", summarise, grp.med=median(V10/1000))

rbreak_table<-rbind(merge(mu,med),c("all breaks",mean(rbreak$V10/1000),median(rbreak$V10/1000)))
rbreak_table$grp.mean<-format(round(as.numeric(rbreak_table$grp.mean), 2), nsmall = 2)
rbreak_table$grp.med<-format(round(as.numeric(rbreak_table$grp.med), 2), nsmall = 2)

stable.p <- ggtexttable(rbreak_table, rows=NULL, cols = c("Group","Group Mean", "Group Median"), 
                        theme = ttheme("mOrange"))+ theme_void()+labs(
                          title = "Recombination Breakpoint Statistics",
                          subtitle = "Resolution of recombination breakpoints are scaled by 1kb. 
Hotspot and non-hotspot sizes were compared in a Mann-Whitney U 
two-sided test (p-value = 0.033)."
                        )

ggplot(rbreak, aes(x=(V10)/1000,fill=significant,color=significant)) +
  geom_histogram(position="identity", alpha=0.5)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=significant),linetype="dashed")+
  labs(y="Number of Recombination Breakpoints", x="Resolution of Recombination Breakpoints (1kb)")
wilcox.test((V10) ~ significant, data=rbreak, alternative="two.sided")
#significant grp.mean
#1     hotspot    6.455
#2 non-hotspot   11.230
#W = 1761750, p-value = 0.03299

mu_log <- ddply(rbreak, "significant", summarise, grp.mean=mean(log10(V10)))
rbreak_grps<-ggplot(rbreak, aes(x=log10(V10),fill=significant,color=significant)) +
  geom_histogram(position="identity", alpha=0.5)+
  geom_vline(data=mu_log, aes(xintercept=grp.mean, color=significant),linetype="dashed")+
  labs(y="Number of Recombination Breakpoints", x="Resolution of Recombination Breakpoints 
       Log10(base pairs)")+
  theme(legend.position = c(0.2, 0.8))

wilcox.test(log10(V10) ~ significant, data=rbreak, alternative="two.sided")
rbreak_all<-ggplot(rbreak, aes(x=log10(V10))) +
  geom_histogram(position="identity", alpha=0.5)+
  labs(y="Number of Recombination Breakpoints", x="Resolution of Recombination Breakpoints 
       Log10(base pairs)")

ggarrange(rbreak_grps, rbreak_all, stable.p, 
          ncol = 1, nrow = 3,  align = "h",
          heights=c(2.2,2.2,2),
          common.legend = FALSE)

########## Plot Recombination Breakpoints: Hotspots vs Non-Hotspots
rbreak_dens <- data.frame(lapply(rbreak_dens, function(x) {gsub("PTT_", "", x)}))
rbreak_dens <- data.frame(lapply(rbreak_dens, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))

rbreak_dens[,2:6] <- data.frame(lapply(rbreak_dens[,2:6], function(x) {as.numeric(as.character(x))}))

rbreak_dens$V1 <- factor(rbreak_dens$V1,levels=c(goodChrOrder,
                                     "0-1_contig_32","0-1_contig_39",
                                     "0-1_contig_40","0-1_contig_43",
                                     "0-1_contig_44","0-1_contig_47",
                                     "0-1_contig_48","0-1_contig_50",
                                     "0-1_contig_51","0-1_contig_52",
                                     "0-1_contig_53","0-1_contig_54",
                                     "0-1_contig_55","0-1_contig_56",
                                     "0-1_contig_57","0-1_contig_59",
                                     "0-1_contig_60","0-1_contig_61",
                                     "0-1_contig_62","0-1_contig_63",
                                     "0-1_contig_65","0-1_contig_67",
                                     "0-1_contig_68","0-1_contig_71",
                                     "0-1_contig_72","0-1_contig_73",
                                     "0-1_contig_74","0-1_contig_75",
                                     "0-1_contig_78","0-1_contig_79",
                                     "0-1_contig_80","0-1_contig_81",
                                     "0-1_contig_82","0-1_contig_84"
))

colnames(rbreak_dens)<-c("chrom","start","stop","breaks","ppois","Bonferroni","significant")

ggplot(rbreak_dens, aes(y=breaks)) +
  geom_linerange(aes(xmin = start, xmax = stop, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("red", "gray"))+
  facet_grid(~chrom, scales = "free_x", space = "free_x")+theme_bw()+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data=rbreak_dens, aes(x=start, y=breaks)) +
  geom_line()+ 
  facet_grid(~chrom, scales = "free_x", space = "free_x")+
  theme(panel.margin = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_linerange(aes(xmin = start, xmax = stop,y=-1, colour=significant),size = 5,show.legend = FALSE)+
  scale_color_manual(values=c("#FF0000FF", "#FFFFFF00"))

rhot<-rbreak_dens[ which(rbreak_dens$significant=='hotspot' ),]
write.table(rhot, file = "all.recombhotspots.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

############################################################################
###################### Genomic Features and Segregation Distortion #########
############################################################################
####format gene density
genedensity<-read.table("PTTxPTM_tracks.genedensity.bed", header = FALSE,sep = "\t")
colnames(genedensity)<-c("chrom","start","stop","genecount")
genedensity<-genedensity[grepl("chrPTT|m86|contig", genedensity$chrom), , drop = FALSE]
genedensity<-data.frame(lapply(genedensity, function(x) {gsub("PTT_", "", x)}))
genedensity<- data.frame(lapply(genedensity, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
genedensity[,2:4] <- data.frame(lapply(genedensity[,2:4], function(x) {as.numeric(as.character(x))}))
genedensity_genome<-genedensity[grepl("chr|mitochondria$", genedensity$chrom), , drop = FALSE]
genedensity_contig<-genedensity[grepl("0-1_contig", genedensity$chrom), , drop = FALSE]
genedensity_genome$chrom <- factor(genedensity_genome$chrom,levels=goodChrOrder)

####format repeats
#repeats<-read.table(PTMxPTT.allrepeat.04292020.bed", header = FALSE,sep = "\t")
#colnames(repeats)<-c("chrom","start","stop","repeat","repeat_type")
#repeats<-repeats[grepl("chrPTT|m86|contig", repeats$chrom), , drop = FALSE]
#repeats<-data.frame(lapply(repeats, function(x) {gsub("PTT_", "", x)}))
#repeats<- data.frame(lapply(repeats, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
#repeats[,2:3] <- data.frame(lapply(repeats[,2:3], function(x) {as.numeric(as.character(x))}))
#repeats_genome<-repeats[grepl("chr|mitochondria$", repeats$chrom), , drop = FALSE]
#repeats_contig<-repeats[grepl("0-1_contig", repeats$chrom), , drop = FALSE]

####format effectors
effectors<-read.table("PTTxPTM_tracks.effectors.bed", header = FALSE,sep = "\t")
colnames(effectors)<-c("chrom","start","stop")
effectors<-effectors[grepl("chrPTT|m86|contig", effectors$chrom), , drop = FALSE]
effectors<-data.frame(lapply(effectors, function(x) {gsub("PTT_", "", x)}))
effectors<- data.frame(lapply(effectors, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
effectors[,2:3] <- data.frame(lapply(effectors[,2:3], function(x) {as.numeric(as.character(x))}))
effectors$type<-"effector"
effectors_genome<-effectors[grepl("chr|mitochondria$", effectors$chrom), , drop = FALSE]
effectors_contig<-effectors[grepl("0-1_contig", effectors$chrom), , drop = FALSE]
effectors_genome$chrom <- factor(effectors_genome$chrom,levels=goodChrOrder)
####format biosynthetic clusters
bsc<-read.table("PTTxPTM_tracks.bsc.bed", header = FALSE,sep = "\t")
colnames(bsc)<-c("chrom","start","stop","region")
bsc<-bsc[grepl("chrPTT|m86|contig", bsc$chrom), , drop = FALSE]
bsc<-data.frame(lapply(bsc, function(x) {gsub("PTT_", "", x)}))
bsc<- data.frame(lapply(bsc, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
bsc[,2:3] <- data.frame(lapply(bsc[,2:3], function(x) {as.numeric(as.character(x))}))
bsc$types<-"bsc"
bsc_genome<-bsc[grepl("chr|mitochondria$", bsc$chrom), , drop = FALSE]
bsc_contig<-bsc[grepl("0-1_contig", bsc$chrom), , drop = FALSE]
bsc_genome$chrom <- factor(bsc_genome$chrom,levels=goodChrOrder)
####format mat loci
mat<-read.table("matloci.bed", header = FALSE,sep = "\t")
colnames(mat)<-c("chrom","start","stop")
mat<-mat[grepl("chrPTT|m86|contig", mat$chrom), , drop = FALSE]
mat<-data.frame(lapply(mat, function(x) {gsub("PTT_", "", x)}))
mat[,2:3] <- data.frame(lapply(mat[,2:3], function(x) {as.numeric(as.character(x))}))
mat$types<-"mat"
mat$chrom <- factor(mat$chrom,levels=goodChrOrder)
####format Presence/Absence Variation
pav<-read.table("PTTvsPTM.nuclearchrmitochondria.shared.txt", header = FALSE,sep = "\t")
colnames(pav)<-c("chrom","start","stop")
pav<-pav[grepl("chrPTT|m86|contig", pav$chrom), , drop = FALSE]
pav<-data.frame(lapply(pav, function(x) {gsub("PTT_", "", x)}))
pav<- data.frame(lapply(pav, function(x) {gsub("0-1_contig_m86","mitochondria", x)}))
pav[,2:3] <- data.frame(lapply(pav[,2:3], function(x) {as.numeric(as.character(x))}))
pav_genome<-pav[grepl("chr|mitochondria$", pav$chrom), , drop = FALSE]
pav_contig<-pav[grepl("0-1_contig", pav$chrom), , drop = FALSE]

pav_genome$types<-"shared"
pav_genome$chrom <- factor(pav_genome$chrom,levels=goodChrOrder)

####format recombination breakpoints
#rbreak_genome<-rbreak_dens[grepl("chr|mito", rbreak_dens$chrom), , drop = FALSE]
#rbreak_contig<-rbreak_dens[grepl("0-1_contig", rbreak_dens$chrom), , drop = FALSE]

####format regions unique to PTTxPTM segregation distortion
#unique<-read.table("uniq.PTTxPTM.bed", header = FALSE,sep = "\t")
#colnames(unique)<-c("chrom","start","stop")
#unique$type<-"unique"
#unique[,2:3] <- data.frame(lapply(unique[,2:3], function(x) {as.numeric(as.character(x))}))
#unique$chrom = factor(unique$chrom,levels=goodChrOrder)
#### Format chromosome names
chrom_names<-c("1","2","3","4","5","6","7","8","9","10","11","12","m")
names(chrom_names)<-goodChrOrder
### Chromosomes and Mitochondria
### Order of Colors:   
p<-ggplot() + 
  geom_hline(yintercept = -0.5,color="white",size = 58)+
  #geom_linerange(data=unique,aes(xmin = start, xmax = stop,y=0,
  #                               colour=type,alpha = 0.5),size=300,show.legend = FALSE)+
  geom_linerange(data=segreg_genome, aes(xmin = start, xmax = stop, y=freq_ptm,
    colour=significant),size = 5,show.legend = FALSE)+
  facet_grid(~chrom,   scales = "free_x",
             labeller = labeller(chrom = chrom_names))+
  scale_color_manual(values=c(
    "bsc"="deepskyblue3",
    "effector"="palegreen3",
    #"hotspot"="red", 
    #"low-simple_repeat"="lightgoldenrodyellow", 
    "mat"="purple",
    #"non-hotspot"="#FFFFFF00",
    "non-significant"="gray",
    "significant"="blue",
    #"TE"="salmon",
    "shared"="darkorchid3"
    #"unique"="light gray"
   ))+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = -0.08,color="gray")+
  theme(panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank())+  
  geom_rect(data=genedensity_genome,aes(xmin = start, xmax = stop,ymin=-0.15,
                                        ymax=-0.1,fill=genecount),show.legend = FALSE)+
  scale_fill_gradient(low="#FFFFFF00",high="gray19")+
  #geom_rect(data=repeats_genome,aes(xmin = start, xmax = stop,ymin=-0.25,
  #                                      ymax=-0.2, colour=repeat_type),show.legend = FALSE)+
  #geom_linerange(data=repeats_genome,aes(xmin = start, xmax = stop,y=-0.23,
  #                                   colour=repeat_type),size=3.5,show.legend = FALSE)+
  geom_rect(data=effectors_genome,aes(xmin = start, xmax = stop,ymin=-0.35,
                                    ymax=-0.3, colour=type),show.legend = FALSE)+
  geom_linerange(data=bsc_genome,aes(xmin = start, xmax = stop,y=-0.43,
                                     colour=types),size=3.5,show.legend = FALSE)+
  geom_linerange(data=mat,aes(xmin = start-10000, xmax = stop+10000,y=-0.5,
                              colour=types),size=5,alpha=1,show.legend = FALSE)+
  #geom_linerange(data=rbreak_genome,aes(xmin = start, xmax = stop,y=-0.6,
  #                                      colour=significant),size=3.5,show.legend = FALSE)+
  geom_linerange(data=pav_genome,aes(xmin = start, xmax = stop,y=-0.7,
                                     colour=types, alpha=0.5),size=3.5,show.legend = FALSE)+
  theme_light()+ theme(panel.margin = unit(0, "lines"),
                    panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
                    #strip.text.x.top = element_text(angle = 90),
                    strip.text.x.top = element_text(size = 20, color = "black",face = "bold"),
                    #axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title=element_text(size=18,face="bold"))+
  scale_y_continuous(limits=c(-0.75,1),name ="            FGOB10Ptm-1 Allele Frequency", 
                   labels=c("","","0","0.5","1.0"))+
  scale_x_continuous(breaks=c(0,1000000,2000000,3000000,4000000,5000000,6000000), 
                     labels=c("0","1","2","3","4","5","6"))+
  xlab("Nuclear and Mitochondrial Genome")

g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
pal <- rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2"))
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  #l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
  #g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- pal[i + 1]
}

pdf("SegregationDistortion_genomicfeatures_nuclmito_PTTxPTM.pdf", width = 21, height=7)
plot(g)
graphics.off()

#####Candidate Regions
candidate_regions<-function(CHROM,range_start,range_stop){
  pal_df<-as.data.frame(rbind(chrom_names,pal))
  ggplot() + 
    geom_hline(yintercept = -0.5,color="white",size = 58)+
    geom_linerange(data=segreg_genome[segreg_genome$chrom==CHROM & 
                                        segreg_genome$start>range_start-1000000 &
                                        segreg_genome$stop<range_stop+1000000,], aes(xmin = start, xmax = stop, y=freq_ptm,
                                                                                     colour=significant),size = 5,show.legend = FALSE)+
    facet_grid(~chrom,   scales = "free_x",
               labeller = labeller(chrom = chrom_names))+
    scale_color_manual(values=c(
      "bsc"="deepskyblue3",
      "effector"="palegreen3",
      "mat"="purple",
      "non-significant"="gray",
      "significant"="blue"
    ))+
    geom_hline(yintercept = 0.5)+
    geom_hline(yintercept = -0.08,color="gray")+
    theme(panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
          strip.text.x.top = element_text(angle = 90),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.y = element_blank())+  
    geom_rect(data=genedensity_genome[genedensity_genome$chrom==CHROM & 
                                        genedensity_genome$start>range_start-1000000 &
                                        genedensity_genome$stop<range_stop+1000000,],aes(xmin = start, xmax = stop,ymin=-0.15,
                                                                                         ymax=-0.1,fill=genecount),show.legend = FALSE)+
    scale_fill_gradient(low="#FFFFFF00",high="gray19")+
    geom_rect(data=effectors_genome[effectors_genome$chrom==CHROM & 
                                      effectors_genome$start>range_start-1000000 &
                                      effectors_genome$stop<range_stop+1000000,],aes(xmin = start, xmax = stop,ymin=-0.35,
                                                                                     ymax=-0.3, colour=type),show.legend = FALSE)+
    geom_linerange(data=bsc_genome[bsc_genome$chrom==CHROM & 
                                     bsc_genome$start>range_start-1000000 &
                                     bsc_genome$stop<range_stop+1000000,],aes(xmin = start, xmax = stop,y=-0.43,
                                                                              colour=types),size=3.5,show.legend = FALSE)+
    geom_linerange(data=pav_genome[pav_genome$chrom==CHROM & 
                                     pav_genome$start>range_start-1000000 &
                                     pav_genome$stop<range_stop+1000000,],aes(xmin = start, xmax = stop,y=-0.7,
                                                                              colour=types, alpha=0.5),size=3.5,show.legend = FALSE)+
    theme_light()+ theme(panel.margin = unit(0, "lines"),
                         panel.border = element_rect(color = "light gray", fill = NA, size = 0.5),
                         strip.text.x.top = element_text(size = 20, color = "black",face = "bold"),
                         strip.background = element_rect(fill=pal_df[2,CHROM]),
                         axis.ticks.x=element_blank(),
                         axis.ticks.y=element_blank(),
                         axis.title=element_text(size=18,face="bold"))+
    scale_y_continuous(limits=c(-0.75,1),name ="            ", 
                       labels=c("","","0","0.5","1.0"))+
    scale_x_continuous(breaks=c(0,1000000,2000000,3000000,4000000,5000000,6000000), 
                       labels=c("0","1","2","3","4","5","6"))+
    xlab(paste(CHROM,":", range_start-1000000, "-", range_stop+1000000))
  
}

chr3_plot<-candidate_regions(CHROM="chr3",range_start=3617625,range_stop=3618276)+
  scale_y_continuous(limits=c(-0.75,1),name ="            FGOB10Ptm-1 Allele Frequency", 
                     labels=c("","","0","0.5","1.0")) ####LD,SD,effectors
chr6_plot<-candidate_regions(CHROM="chr6",range_start=1503851,range_stop=1504863) ####LD,SD,Dxy,Fst,effectors
chr10_plot<-candidate_regions(CHROM="chr10",range_start=1206644,range_stop=1264872) ### SD,bsc,(PTT)

candidate_region_plots<-ggarrange(chr3_plot,chr6_plot,chr10_plot, 
                                  ncol = 3, labels = c("B", "C","D"),
                                  font.label = list(size = 30, family = "Times"))+
  theme(plot.margin = margin(1,0,0,0, "cm"))




pdf("Figure3_SegregationDistortion_genomicfeatures_nuclmito_PTTxPTM.pdf", width = 21, height=14)
print(ggarrange(g, candidate_region_plots, 
          nrow = 2, labels = "A",
          font.label = list(size = 30, family = "Times")))
graphics.off()
