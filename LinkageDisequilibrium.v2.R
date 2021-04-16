#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
#Rscript --vanilla LinkageDisequilibrium.R all.r.ptm.matrix.ld all.r2.ptm.matrix.ld all.r2.ptm.ld all.r.ptm.inter.ld

###install.packages("corrplot")
###library("corrplot")
###library("ggplot2")

library(devtools)
#install_github("jokergoo/ComplexHeatmap", force = TRUE)
library(ComplexHeatmap)
library(circlize)
library(data.table)

r_matrix<-as.matrix(read.table(args[1], header = FALSE,sep = "\t",na.strings = "nan"))
r2_matrix<-as.matrix(read.table(args[2], header = FALSE,sep = "\t",na.strings = "nan"))
r<-as.data.frame(read.table(args[3], header = TRUE))
r_inter<-as.data.frame(read.table(args[4], header = TRUE))


goodChrOrder <- c(1:12,"m","0-1_contig_32","0-1_contig_39",
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
                                     "0-1_contig_82","0-1_contig_84")

#r_matrix<-as.matrix(read.table("all.r.ptm.matrix.ld", header = FALSE,sep = "\t",na.strings = "nan"))
#r2_matrix<-as.matrix(read.table("all.r2.ptm.matrix.ld", header = FALSE,sep = "\t",na.strings = "nan"))
#r<-as.data.frame(read.table("all.r2.ptm.ld", header = TRUE))
#r_inter<-as.data.frame(read.table("all.r.ptm.inter.ld", header = TRUE))
odds_raw<-as.data.frame(read.table("all.inter.odds",sep = " ", header = FALSE))

#######make odds dataframe symmetric##############
odds_rev<-NULL
odds_rev$V1<-odds_raw$V2 
odds_rev$V2<-odds_raw$V1 
odds_rev$V3<-odds_raw$V4 
odds_rev$V4<-odds_raw$V3 
odds_rev$V5<-odds_raw$V8 
odds_rev$V6<-odds_raw$V7 
odds_rev$V7<-odds_raw$V6 
odds_rev$V8<-odds_raw$V5 
odds_rev$V9<-odds_raw$V9 
odds_rev$V10<-odds_raw$V10

odds<-data.frame(rbind(odds_raw,data.frame(odds_rev,stringsAsFactors = FALSE)),stringsAsFactors = FALSE)
odds<-unique(odds[c("V1","V2","V3","V4","V9")])


############################################
############Plot R/R2 values###################
############################################
markers<-read.table("all.markers.bed", header = FALSE,sep = "\t")
markers$V1 <- as.character(markers$V1)
markers$V1[markers$V1 == "0-1_contig_m86"] <- "m"
markers$V1<-gsub("chrPTT_","",markers$V1)
markers$V1<-factor(markers$V1,levels=goodChrOrder, ordered=TRUE)
markers$V2<-factor(markers$V2,levels=c(as.character(markers$V2)) ,ordered=TRUE)
markers_odds_V1<-unique(sort(odds$V3))
markers_odds_V2<-unique(sort(odds$V4))

colnames(r_matrix)<-markers[,2]
rownames(r_matrix)<-markers[,2]
colnames(r2_matrix)<-markers[,2]
rownames(r2_matrix)<-markers[,2]

### Odds Ratio: construct 0 matrix of correct dimensions with row and column names
myMat_odds <- matrix(NA, length(markers_odds_V1), length(markers_odds_V2), 
	dimnames = list(markers_odds_V1, markers_odds_V2))
myMat_odds[as.matrix(odds[c("V3", "V4")])] <- odds[["V9"]]

##### sort chromosomes in marker dataframe and r/r2_matrix 
####markers<-markers[order(markers$V1),]
markers<-markers[order(match(markers[[1]], goodChrOrder)), ]

r_matrix<-r_matrix[,match(markers$V2, colnames(r_matrix))]
r_matrix<-r_matrix[,match(markers$V2, rownames(r_matrix))]
r2_matrix<-r2_matrix[,match(markers$V2, colnames(r2_matrix))]
r2_matrix<-r2_matrix[,match(markers$V2, rownames(r2_matrix))]

markers_odds_V1<-na.omit(markers_odds_V1[match(markers$V2, markers_odds_V1)])
markers_odds_V2<-na.omit(markers_odds_V2[match(markers$V2, markers_odds_V2)])
myMat_odds<-myMat_odds[,match(markers_odds_V2, colnames(myMat_odds))]
myMat_odds<-myMat_odds[match(markers_odds_V1, rownames(myMat_odds)),]

##### use only chromosomes and leave out contigs (don't know where contigs are in relation to chromosomes--inter/intra chromosomal) ####
colsub<-colnames(r_matrix)[!grepl("contig", colnames(r_matrix))]
rowsub<-colnames(r_matrix)[!grepl("contig", rownames(r_matrix))]
colsub_odds<-colnames(myMat_odds)[!grepl("contig", colnames(myMat_odds))]
rowsub_odds<-rownames(myMat_odds)[!grepl("contig", rownames(myMat_odds))]

markersub<-markers[match(rowsub,markers$V2),]
markersub$V1<-factor(markersub$V1,levels=goodChrOrder, ordered=TRUE)
r_matrix<-r_matrix[rowsub,colsub]
r2_matrix<-r2_matrix[rowsub,colsub]
myMat_odds<-myMat_odds[rowsub_odds,colsub_odds]
rowsub_odd<-as.data.frame(cbind(rowsub_odds,rowsub_odds))
rowsub_odds_split <- gsub("(:.*|chrPTT_)", "", rowsub_odds)
colsub_odds_split <- gsub("(:.*|chrPTT_)", "", colsub_odds)

colors<-c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2")

dim(markers)
dim(r_matrix)

####### Plot R ##########
plot.r<-Heatmap(r_matrix,row_split = markersub$V1, column_split=markersub$V1,show_row_names = FALSE,show_column_names = FALSE,
	row_title = NULL, column_title=NULL, 
	heatmap_legend_param = list(title = "R"),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(markersub$V1)),
				       labels_gp = gpar(col = "white", fontsize = 20))),
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(markersub$V1)),
                                       labels_gp = gpar(col = "white", fontsize = 20))),
        na_col = "gray", cluster_rows=FALSE, cluster_columns=FALSE,raster_device = "CairoPNG")

png(filename = "all.r.png",
    width = 1080, height = 1080, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
plot.r
dev.off()

#### Plot R2 ####
plot.r2<-Heatmap(r2_matrix,row_split = markersub$V1, column_split=markersub$V1,show_row_names = FALSE,show_column_names = FALSE,
        row_title = NULL, column_title=NULL,
        heatmap_legend_param = list(title = "R"),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(markersub$V1)),
                                       labels_gp = gpar(col = "white", fontsize = 20))),
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(markersub$V1)),
                                       labels_gp = gpar(col = "white", fontsize = 20))),
        na_col = "gray", cluster_rows=FALSE, cluster_columns=FALSE,raster_device = "CairoPNG")

png(filename = "all.r2.png",
    width = 1080, height = 1080, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
plot.r2
dev.off()

######## Plot Odds Ratio ############
plot.odds<-Heatmap(myMat_odds,row_split = rowsub_odds_split, column_split=colsub_odds_split,show_row_names = FALSE,show_column_names = FALSE,
        row_title = NULL, column_title=NULL,
	col=colorRamp2(c(0, 0.8, 1,3,10,max(odds$V9)), c("purple4","blue", "green4","yellow","orange", "red")),
        heatmap_legend_param = list(title = "R"),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(rowsub_odds_split)),
                                       labels_gp = gpar(col = "white", fontsize = 20))),
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = colors),
                                       labels = c(unique(colsub_odds_split)),
                                       labels_gp = gpar(col = "white", fontsize = 20))),
        na_col = "gray", cluster_rows=FALSE, cluster_columns=FALSE,raster_device = "CairoPNG")

png(filename = "all.odds.png",
    width = 1080, height = 1080, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
plot.odds
dev.off()

####### Plot R/R2 histogram #########
##### R interchromosomal
r_inter<-r_inter[r_inter$SNP_A %in% r_inter$SNP_A[!grepl("contig", r_inter$SNP_A)],]
r_inter<-r_inter[r_inter$SNP_B %in% r_inter$SNP_B[!grepl("contig", r_inter$SNP_B)],]
r_inter$R2<-r_inter$R^2

png(filename = "all.inter.r.hist.png",
    width = 800, height = 800, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
hist(r_inter$R)
dev.off()


##### R2 interchromosomal 
png(filename = "all.inter.r2.hist.png",
    width = 800, height = 800, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
hist(log10(r_inter$R2))
dev.off()

###### Odds Ratio interchromosomal
png(filename = "all.inter.odds.hist.png",
    width = 800, height = 800, units = "px", pointsize = 12,
    bg = "white",  res = NA,type = c("cairo"))
#hist(na.omit(odds_raw$V9))
hist(log10(na.omit(odds_raw$V9)))
dev.off()

########### Top 5% of R/R2/Odds Ratio values ##############
odds_top.1perc<-odds_raw[odds_raw$V9 >quantile(odds_raw$V9,prob=0.999),]
min(odds_top.1perc$V9) #30.1652
write.table(odds_top.1perc, "odds_top.1perc.bed", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE)
