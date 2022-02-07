#module load R/3.5.3
library(circlize) 
library(migest)
library(dplyr)
library(ComplexHeatmap)

###Files
PTTxPTM<-read.table("PTTxPTM_tracks.nuclmito.txt", header=FALSE, sep="\t")
PTTxPTM_chords<-read.table("PTTvsPTM.nuclearchrmitochondria.masked.show-coords.txt", header=FALSE, sep="\t")
PTTxPTM_inversions_chords<-read.table("PTTvsPTM.nuclearchr_mito.masked.inversions.show-coords", header=FALSE, sep="\t")
PTTxPTM_r2_chords<-read.table("r_inter.1perc.trim.rpos.genic.show-coords.txt", header=FALSE, sep="\t")
PTTxPTM_odds_chords<-read.table("odds_top.1perc.trim.ptmptt.genic_nogene.show-coords.txt", header=FALSE, sep="\t")

PTTxPTM_inversions<-read.table("PTTvsPTM.nuclear_mito.masked.inversions.tracks.txt", header=FALSE, sep="\t")
PTTxPTM_genes<-read.table("PTTxPTM_tracks.nuclmito_genedensity.bed", header=FALSE, sep="\t")
PTTxPTM_repeats<-read.table("PTTxPTM_tracks.nuclmito_repeat.bed", header=FALSE, sep="\t")
PTTxPTM_effectors<-read.table("PTTxPTM_tracks.nuclmito_effectors.bed", header=FALSE, sep="\t")
PTTxPTM_bsc<-read.table("PTTxPTM_tracks.nuclmito_bsc.bed", header=FALSE, sep="\t")
PTTxPTM_mat<-read.table("matloci.bed", header=FALSE, sep="\t")
PTTxPTM_rhot<-read.table("all.recombhotspots.txt", header=FALSE, sep="\t")
PTTxPTM_shared<-read.table("PTTvsPTM.nuclearchrmitochondria.shared.txt", header=FALSE, sep="\t")


###recombination breakpoints (red)
###segregation distortion (blue/red histogram)


colors<-c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2",
          "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2")
PTTxPTM$colors<-rev(colors)
colnames(PTTxPTM)<-c("chr", "start", "stop","colors")
goodPTTOrder <- c(paste("chrPTT_",c(1:12),sep=""),"chrPTT_m")
goodPTMOrder <- c("chrPTM_m",paste("chrPTM_",c(12:1),sep=""))

genome<-"PTT"
exclude_genome<-"PTM"
goodGenomeOrder<-goodPTTOrder

##cdm_res =chordDiagram(PTTxPTM, annotationTrack = "grid", order=c(goodPTMOrder, goodPTTOrder))
PTTxPTM$chr<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM$chr)
PTTxPTM$chr<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM$chr)
PTTxPTM$chr <- factor(PTTxPTM$chr, levels = c(goodPTMOrder, goodPTTOrder))
PTTxPTM<-PTTxPTM[match(c(goodPTMOrder, goodPTTOrder), PTTxPTM$chr),]
PTTxPTM<-PTTxPTM[!grepl(exclude_genome, PTTxPTM$chr),]

PTTxPTM_odds_chords<-PTTxPTM_odds_chords[!grepl(exclude_genome, PTTxPTM_odds_chords$V1),]

PTTxPTM_chords$V2<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_chords$V2)
PTTxPTM_chords$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_chords$V1)

PTTxPTM_inversions$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_inversions$V1)
PTTxPTM_inversions$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_inversions$V1)
PTTxPTM_inversions<-PTTxPTM_inversions[!grepl(exclude_genome, PTTxPTM_inversions$V1),]

PTTxPTM_genes$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_genes$V1)
PTTxPTM_genes$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_genes$V1)
PTTxPTM_genes<-PTTxPTM_genes[!grepl(exclude_genome, PTTxPTM_genes$V1),]

PTTxPTM_repeats$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_repeats$V1)
PTTxPTM_repeats$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_repeats$V1)
PTTxPTM_repeats<-PTTxPTM_repeats[!grepl(exclude_genome, PTTxPTM_repeats$V1),]

PTTxPTM_effectors$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_effectors$V1)
PTTxPTM_effectors$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_effectors$V1)
PTTxPTM_effectors<-PTTxPTM_effectors[!grepl(exclude_genome, PTTxPTM_effectors$V1),]

PTTxPTM_bsc<-PTTxPTM_bsc[!grepl(exclude_genome, PTTxPTM_bsc$V1),]
PTTxPTM_mat<-PTTxPTM_mat[!grepl(exclude_genome, PTTxPTM_mat$V1),]

PTTxPTM_rhot$V1<-gsub("chr","chrPTT_",PTTxPTM_rhot$V1)
PTTxPTM_rhot$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_rhot$V1)
PTTxPTM_rhot$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_rhot$V1)
PTTxPTM_rhot<-PTTxPTM_rhot[!grepl(exclude_genome, PTTxPTM_rhot$V1),]
PTTxPTM_rhot<-PTTxPTM_rhot[!grepl("contig", PTTxPTM_rhot$V1),]

PTTxPTM_shared$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_shared$V1)
PTTxPTM_shared$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_shared$V1)
PTTxPTM_shared<-PTTxPTM_shared[!grepl(exclude_genome, PTTxPTM_shared$V1),]
##############################################
########### Genome Compare Plot #############
##############################################

#png("Circlize_plot_may282020.png", type="cairo", height=1000,width=2000,units = "px", pointsize = 12,)
#pdf("Circlize_plot_may282020.pdf", height=16, width=23)
pdf(paste0("Circlize_plot_oddsratio_",genome,"july282020.pdf"), height=16, width=23)

#pdf("Aquavit.circlize_plot_r2.pdf", height=16, width=23)
par(mar=rep(0,4))
circos.clear()
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.05), start.degree = 90, gap.degree =.5,
          canvas.xlim=c(-2.2, 2.2),   # bigger canvas?
           canvas.ylim=c(-1.2, 1.2))   # bigger canvas?
### Sector details
circos.initialize(factors = PTTxPTM$chr, xlim = cbind(PTTxPTM$start, PTTxPTM$stop))
### Plot sectors
PTTPTM_res =circos.trackPlotRegion(ylim = c(-1, 0), factors = goodGenomeOrder, track.height=0.2,
   #panel.fun for each sector
   panel.fun = function(x, y) {
   #select details of current sector
   name = gsub(paste0("chr",genome,"_"),"",get.cell.meta.data("sector.index"))
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
   ylim = get.cell.meta.data("ylim")
   #plot main sector
   circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2],col = PTTxPTM$colors[i])
   #plot country labels
   circos.text(x=mean(xlim), y=-0.5, labels=name, cex = 2,col = "black", facing = "inside",niceFacing = TRUE)
                       })

######### Presence/Absence Variation ##########
for(i in seq_len(nrow(PTTxPTM_shared))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_shared[i,"V2"], ybottom=ylim[1]+3.5, xright=PTTxPTM_shared[i,"V3"], ytop=ylim[1]+4,
               sector.index=PTTxPTM_shared$V1[i],
               col = "darkorchid3", border=NA)
}
######### Gene Density ##########
gene_col=colorRamp2(c(0,100,200,300),
                    c("white","gray","dimgray","black"), transparency = 0.5)
for(i in seq_len(nrow(PTTxPTM_genes))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_genes[i,"V2"], ybottom=ylim[1]+1.1, xright=PTTxPTM_genes[i,"V3"], ytop=ylim[1]+1.5, 
               sector.index=PTTxPTM_genes$V1[i],
               col = gene_col(PTTxPTM_genes$V4[i]), border=NA)
}
########## Repeats #############
PTTxPTM_repeats$colors<-PTTxPTM_repeats$V5
PTTxPTM_repeats$colors<-gsub("TE","salmon",PTTxPTM_repeats$colors)
PTTxPTM_repeats$colors<-gsub("low-simple_repeat","lightgoldenrodyellow",PTTxPTM_repeats$colors)
for(i in seq_len(nrow(PTTxPTM_repeats))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_repeats[i,"V2"], ybottom=ylim[1]+1.5, xright=PTTxPTM_repeats[i,"V3"], ytop=ylim[1]+2, 
               sector.index=PTTxPTM_repeats$V1[i],
               col = PTTxPTM_repeats$colors[i], border=NA)
}

############# Effectors ###########
for(i in seq_len(nrow(PTTxPTM_effectors))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_effectors[i,"V2"], ybottom=ylim[1]+2.1, xright=PTTxPTM_effectors[i,"V3"], ytop=ylim[1]+2.5, 
               sector.index=PTTxPTM_effectors$V1[i],
               col = "palegreen3", border = "palegreen3")
}

########### Biosynthetic Clusters ##########
for(i in seq_len(nrow(PTTxPTM_bsc))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_bsc[i,"V2"], ybottom=ylim[1]+2.5, xright=PTTxPTM_bsc[i,"V3"], ytop=ylim[1]+3, 
               sector.index=PTTxPTM_bsc$V1[i],
               col = "deepskyblue3", border = "deepskyblue3")
}

############# Inversions ###########
#for(i in seq_len(nrow(PTTxPTM_inversions))) {
#   ylim = get.cell.meta.data("ylim")
#   circos.rect(xleft=PTTxPTM_inversions[i,"V2"], ybottom=ylim[1]+3, xright=PTTxPTM_inversions[i,"V3"], ytop=ylim[1]+3.5,
#               sector.index=PTTxPTM_inversions$V1[i],
#               col = "royalblue", border = "royalblue")
#}

############# Mat Loci ###########
for(i in seq_len(nrow(PTTxPTM_mat))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_mat[i,"V2"], ybottom=ylim[1]+3.1, xright=PTTxPTM_mat[i,"V3"], ytop=ylim[1]+3.5,
               sector.index=PTTxPTM_mat$V1[i],
               col = "darkslateblue", border = "darkslateblue")
}

############# Recombination Hotspots ###########
#for(i in seq_len(nrow(PTTxPTM_rhot))) {
#   ylim = get.cell.meta.data("ylim")
#   circos.rect(xleft=PTTxPTM_rhot[i,"V2"], ybottom=ylim[1]+4, xright=PTTxPTM_rhot[i,"V3"], ytop=ylim[1]+4.5,
#               sector.index=PTTxPTM_rhot$V1[i],
#               col = "red", border = "red")
#}

### Plot links #########
#for(k in 1:nrow(PTTxPTM_chords)){
   ###i,j reference of flow matrix
#   i<-match(PTTxPTM_chords$V1[k],PTTxPTM$chr)
#   j<-match(PTTxPTM_chords$V2[k],PTTxPTM$chr)
   
   ###plot link
#   circos.link(sector.index1=PTTxPTM$chr[i], point1=c(PTTxPTM_chords$V3[k], PTTxPTM_chords$V4[k]),
#               sector.index2=PTTxPTM$chr[j], point2=c(PTTxPTM_chords$V5[k], PTTxPTM_chords$V6[k]),
#               col = add_transparency("grey", transparency = 0.9))
   
   ###update sum1 and sum2 for use when plotting the next link
#   PTTxPTM_chords$V3[k] = PTTxPTM_chords$V3[k] 
#   PTTxPTM_chords$V4[k] = PTTxPTM_chords$V4[k]
#   PTTxPTM_chords$V5[k] = PTTxPTM_chords$V5[k] 
#   PTTxPTM_chords$V6[k] = PTTxPTM_chords$V6[k]
#}

#### Plot R2 ###########
### Only Plot potential DMIs
#PTTxPTM_r2_chords<-subset(PTTxPTM_r2_chords,V7<0)

#for(k in 1:nrow(PTTxPTM_r2_chords)){
   #i,j reference of flow matrix
#   i<-match(PTTxPTM_r2_chords$V1[k],PTTxPTM$chr)
#   j<-match(PTTxPTM_r2_chords$V2[k],PTTxPTM$chr)

   #plot link
#   circos.link(sector.index1=PTTxPTM$chr[i], point1=c(PTTxPTM_r2_chords$V3[k], PTTxPTM_r2_chords$V4[k]),
#               sector.index2=PTTxPTM$chr[j], point2=c(PTTxPTM_r2_chords$V5[k], PTTxPTM_r2_chords$V6[k]),
#               col = add_transparency(PTTxPTM$colors[i], transparency = 0.4))

   #update sum1 and sum2 for use when plotting the next link
#   PTTxPTM_r2_chords$V3[k] = PTTxPTM_r2_chords$V3[k]
#   PTTxPTM_r2_chords$V4[k] = PTTxPTM_r2_chords$V4[k]
#   PTTxPTM_r2_chords$V5[k] = PTTxPTM_r2_chords$V5[k]
#   PTTxPTM_r2_chords$V6[k] = PTTxPTM_r2_chords$V6[k]
#}

#### Plot Odds Ratio ################
for(k in 1:nrow(PTTxPTM_odds_chords)){
   ###i,j reference of flow matrix
   i<-match(PTTxPTM_odds_chords$V1[k],PTTxPTM$chr)
   j<-match(PTTxPTM_odds_chords$V2[k],PTTxPTM$chr)

   ###plot link
   circos.link(sector.index1=PTTxPTM$chr[i], point1=c(PTTxPTM_odds_chords$V3[k], PTTxPTM_odds_chords$V4[k]),
               sector.index2=PTTxPTM$chr[j], point2=c(PTTxPTM_odds_chords$V5[k], PTTxPTM_odds_chords$V6[k]),
               col = add_transparency(PTTxPTM$colors[i], transparency = 0.4) )

   ###update sum1 and sum2 for use when plotting the next link
   PTTxPTM_odds_chords$V3[k] = PTTxPTM_odds_chords$V3[k]
   PTTxPTM_odds_chords$V4[k] = PTTxPTM_odds_chords$V4[k]
   PTTxPTM_odds_chords$V5[k] = PTTxPTM_odds_chords$V5[k]
   PTTxPTM_odds_chords$V6[k] = PTTxPTM_odds_chords$V6[k]
}




############ PTT vs PTM genomes ##############
#highlight.sector(goodPTTOrder,  track.index = 1, col  = "firebrick3",padding = c(0.1, 0, -2, 0),lwd=3,
#                 text = "PTT", cex = 0.8, text.col = "white", niceFacing = TRUE)

#highlight.sector(goodPTMOrder, track.index = 1, col  = "mediumblue",padding = c(0.1, 0, -2, 0), lwd=3,
#                 text = "PTM", cex = 0.8, text.col = "white", niceFacing = TRUE)


############# Legend ######################
lgd_gene = Legend(at = c(0,100,200,300), col_fun = gene_col, 
                   title_position = "topleft", title = "gene density (50kb window)")

lgd_repeat = Legend(at = c("transposable element", "low-simple_repeat"), 
                    legend_gp = gpar(fill = c("salmon","lightgoldenrodyellow")),
                    title_position = "topleft", title = "repeat class")
lgd_virulence = Legend(at=c("effector", "biosynthetic cluster"), 
                    legend_gp = gpar(fill = c("palegreen3", "deepskyblue3")),
                    title_position = "topleft", title = "putative pathogenicity factors")

#lgd_inversion = Legend(at = "inversion",
#                    legend_gp = gpar(fill = "royalblue"),
#                    title_position = "topleft", title = "Genome Rearrangement")
lgd_mat = Legend(at = "mat loci",
                    legend_gp = gpar(fill = "darkslateblue"),
                    title_position = "topleft", title = "mating genes")

#lgd_rhot = Legend(at = "hotspots",
#                    legend_gp = gpar(fill = "red"),
#                    title_position = "topleft", title = "recombination")
lgd_pav = Legend(at = c("shared","non-shared"),
                    legend_gp = gpar(fill = c("darkorchid3","white")),
                    border=c("darkorchid3","gray85"),
                    title_position = "topleft", title = "synteny")

lgd_list_vertical = packLegend(lgd_gene,lgd_repeat,lgd_virulence,lgd_mat,lgd_pav)

draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))


dev.off()
