library(circlize) 
library(migest)
library(dplyr)
library(ComplexHeatmap)

###Files
PTTxPTM<-read.table("PTTxPTM_tracks.nuclmito.txt", header=FALSE, sep="\t")
PTTxPTM_chords<-read.table("PTTvsPTM.nuclearchrmitochondria.masked.show-coords.orthologous.txt", header=FALSE, sep="\t")

PTTxPTM_genes<-read.table("PTTvsPTM_tracks.genes.orthologous.bed", header=FALSE, sep="\t")
PTTxPTM_effectors<-read.table("PTTvsPTM_tracks.effector.orthologous.bed", header=FALSE, sep="\t")
PTTxPTM_bsc<-read.table("PTTvsPTM_tracks.bsc.orthologous.bed", header=FALSE, sep="\t")
PTTxPTM_mat<-read.table("matloci.bed", header=FALSE, sep="\t")
PTTxPTM_pav<-read.table("PTTxPTM_tracks.nuclmito_pav.bed", header=FALSE, sep="\t")
PTTxPTM_shared<-read.table("PTTvsPTM.nuclearchrmitochondria.shared.txt", header=FALSE, sep="\t")
###segregation distortion (blue/red histogram)


colors<-c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2",
          "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#92C5DE","#3288BD", "#2166AC", "#5E4FA2")
PTTxPTM$colors<-rev(colors)
colnames(PTTxPTM)<-c("chr", "start", "stop","colors")
goodPTTOrder <- c(paste("chrPTT_",c(1:12),sep=""),"chrPTT_m")
goodPTMOrder <- c("chrPTM_m",paste("chrPTM_",c(12:1),sep=""))
#cdm_res =chordDiagram(PTTxPTM, annotationTrack = "grid", order=c(goodPTMOrder, goodPTTOrder))
PTTxPTM$chr<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM$chr)
PTTxPTM$chr<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM$chr)
PTTxPTM$chr <- factor(PTTxPTM$chr, levels = c(goodPTMOrder, goodPTTOrder))
PTTxPTM<-PTTxPTM[match(c(goodPTMOrder, goodPTTOrder), PTTxPTM$chr),]

PTTxPTM_chords$V2<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_chords$V2)
PTTxPTM_chords$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_chords$V1)


PTTxPTM_genes$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_genes$V1)
PTTxPTM_genes$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_genes$V1)


PTTxPTM_effectors$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_effectors$V1)
PTTxPTM_effectors$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_effectors$V1)

PTTxPTM_pav$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_pav$V1)
PTTxPTM_pav$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_pav$V1)

PTTxPTM_shared$V1<-gsub("mitochondrial_contig","chrPTM_m",PTTxPTM_shared$V1)
PTTxPTM_shared$V1<-gsub("0-1_contig_m86","chrPTT_m",PTTxPTM_shared$V1)
##############################################
########### Genome Compare Plot #############
##############################################

pdf("Circlize_plot_Feb072022.pdf", height=16, width=23)
par(mar=rep(0,4))
circos.clear()
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.05), start.degree = 90, gap.degree =.5,
          canvas.xlim=c(-2.2, 2.2),   # bigger canvas?
           canvas.ylim=c(-1.2, 1.2))   # bigger canvas?
### Sector details
circos.initialize(factors = PTTxPTM$chr, xlim = cbind(PTTxPTM$start, PTTxPTM$stop))
### Plot sectors
PTTPTM_res =circos.trackPlotRegion(ylim = c(-1, 0), factors = PTTxPTM$chr, track.height=0.07,
   #panel.fun for each sector
   panel.fun = function(x, y) {
   #select details of current sector
   name = gsub("chrPT._","",get.cell.meta.data("sector.index"))
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
   ylim = get.cell.meta.data("ylim")
   #plot main sector
   circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2],col = PTTxPTM$colors[i])
   #plot country labels
   circos.text(x=mean(xlim), y=-0.5, labels=name, cex = 0.6,col = "black", facing = "inside",niceFacing = TRUE)
                       })


for(i in seq_len(nrow(PTTxPTM_shared))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_shared[i,"V2"], ybottom=ylim[1]+3.5, xright=PTTxPTM_shared[i,"V3"], ytop=ylim[1]+4,
               sector.index=PTTxPTM_shared$V1[i],
               col = "darkorchid3", border=NA)
}

######### Gene Density ##########
#gene_col=colorRamp2(c(0,100,200,300),
#                    c("white","gray","dimgray","black"), transparency = 0.5)
#for(i in seq_len(nrow(PTTxPTM_genes))) {
#   ylim = get.cell.meta.data("ylim")
#   circos.rect(xleft=PTTxPTM_genes[i,"V2"], ybottom=ylim[1]+1.1, xright=PTTxPTM_genes[i,"V3"], ytop=ylim[1]+1.5, 
#               sector.index=PTTxPTM_genes$V1[i],
#               col = gene_col(PTTxPTM_genes$V4[i]), border=NA)
#}

for(i in seq_len(nrow(PTTxPTM_genes))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_genes[i,"V2"], ybottom=ylim[1]+1.1, xright=PTTxPTM_genes[i,"V3"], ytop=ylim[1]+1.5,
               sector.index=PTTxPTM_genes$V1[i],
               col = "dimgray", border=NA)
}


############# Effectors ###########
for(i in seq_len(nrow(PTTxPTM_effectors))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_effectors[i,"V2"], ybottom=ylim[1]+2, xright=PTTxPTM_effectors[i,"V3"], ytop=ylim[1]+2.5, 
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


############# Mat Loci ###########
for(i in seq_len(nrow(PTTxPTM_mat))) {
   ylim = get.cell.meta.data("ylim")
   circos.rect(xleft=PTTxPTM_mat[i,"V2"], ybottom=ylim[1]+3, xright=PTTxPTM_mat[i,"V3"], ytop=ylim[1]+3.5,
               sector.index=PTTxPTM_mat$V1[i],
               col = "darkslateblue", border = "darkslateblue")
}

### Plot links #########
for(k in 1:nrow(PTTxPTM_chords)){
   ###i,j reference of flow matrix
   i<-match(PTTxPTM_chords$V1[k],PTTxPTM$chr)
   j<-match(PTTxPTM_chords$V2[k],PTTxPTM$chr)
   
   ###plot link
   circos.link(sector.index1=PTTxPTM$chr[i], point1=c(PTTxPTM_chords$V3[k], PTTxPTM_chords$V4[k]),
               sector.index2=PTTxPTM$chr[j], point2=c(PTTxPTM_chords$V5[k], PTTxPTM_chords$V6[k]),
               col = add_transparency(PTTxPTM$colors[i], transparency = 0.9))
   
   ###update sum1 and sum2 for use when plotting the next link
   PTTxPTM_chords$V3[k] = PTTxPTM_chords$V3[k] 
   PTTxPTM_chords$V4[k] = PTTxPTM_chords$V4[k]
   PTTxPTM_chords$V5[k] = PTTxPTM_chords$V5[k] 
   PTTxPTM_chords$V6[k] = PTTxPTM_chords$V6[k]
}


############ PTT vs PTM genomes ##############
highlight.sector(goodPTTOrder,  track.index = 1, col  = "firebrick3",padding = c(0.1, 0, -2, 0),lwd=3,
                 text = "PTT", cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(goodPTMOrder, track.index = 1, col  = "mediumblue",padding = c(0.1, 0, -2, 0), lwd=3,
                 text = "PTM", cex = 0.8, text.col = "white", niceFacing = TRUE)


############# Legend ######################
lgd_gene = Legend(at = "gene",
                   legend_gp = gpar(fill = "dimgray"),
                   title_position = "topleft", title = "gene")

lgd_virulence = Legend(at=c("effector", "biosynthetic cluster"), 
                    legend_gp = gpar(fill = c("palegreen3", "deepskyblue3")),
                    title_position = "topleft", title = "putative virulence factors")

lgd_mat = Legend(at = "mat loci",
                    legend_gp = gpar(fill = "darkslateblue"),
                    title_position = "topleft", title = "mating genes")

lgd_pav = Legend(at = c("shared","non-shared"),
                    legend_gp = gpar(fill = c("darkorchid3","white")),
		    border=c("darkorchid3","gray85"),
                    title_position = "topleft", title = "synteny")

lgd_list_vertical = packLegend(lgd_gene,lgd_virulence,lgd_mat,lgd_pav)

draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))


graphics.off()


#qpdf::pdf_combine(input = c("Figure1A_lesiontypes.pdf", 
#                            "Figure1B_Circlize_plot_Feb072022.pdf"),
#                  output = "output.pdf")
