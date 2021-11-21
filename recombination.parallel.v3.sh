#!/bin/bash

#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=1PTTxPTMrecomb
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=1G
#  write std out and std error to these files
#SBATCH --error=PTTxPTMrecomb.%J.err
#SBATCH --output=PTTxPTMrecomb.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=global

#################################################################################################################
############################Create a r2 and r matrix#############################################################
#################################### Plink ######################################################################
#grep "##" ../32.PTMPTT.flt.vcf >header
#awk -v OFS="\t" '{print $1,$2,$1":"$2"-"$3,"T","A",".",".",".","GT:AD:DP:GQ:PL"}' all.recomblocks.bed >ptt.columns
#sed 's/PTT_/ptt_/g' ptt.columns| sed 's/0-1_contig/ptm/g'|tail -n +2 >ptm.columns
#cat ptt.columns ptm.columns \
#|sed  "1s/.*/#CHROM POS ID	 REF ALT QUAL FILTER INFO FORMAT/" \
#|awk '{$1=$1}1' OFS="\t" >columns
#sed 's/PTT/0\/0/g' all.recomblocks.vcf.tmp |sed 's/PTM/0\/1/g' | sed 's/\./\.\/\./g'|paste columns - | cat header - >all.recomblocks.ptm.vcf
#sed 's/PTT/0\/1/g' all.recomblocks.vcf.tmp |sed 's/PTM/0\/0/g' | sed 's/\./\.\/\./g'|paste columns - | cat header - >all.recomblocks.ptt.vcf

#vcftools --vcf all.recomblocks.ptm.vcf --max-missing-count 100 --out all.recomblocks.ptm --recode
#vcftools --vcf all.recomblocks.ptt.vcf --max-missing-count 100 --out all.recomblocks.ptt --recode
#####check if error such that PTT and PTM alles in an isolate overlap
#grep 'chrPTT_.*chrPTT_.*PTM.*chrPTT_.*PTT' *.ranges.bed
#grep 'contig.*contig.*PTM.*contig.*PTT' *.ranges.bed

#Create a r2 and r matrix
##Plink
#/data/biosoftware/plink/v1.90beta4/plink --r square --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r.ptm.matrix --allow-extra-chr
#/data/biosoftware/plink/v1.90beta4/plink --r2 square --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r2.ptm.matrix --allow-extra-chr

#/data/biosoftware/plink/v1.90beta4/plink --r square --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r.ptt.matrix --allow-extra-chr
#/data/biosoftware/plink/v1.90beta4/plink --r2 square --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r2.ptt.matrix --allow-extra-chr

#/data/biosoftware/plink/v1.90beta4/plink --r2 inter-chr --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r2.ptm --allow-extra-chr
#/data/biosoftware/plink/v1.90beta4/plink --r inter-chr --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r.ptm --allow-extra-chr

#/data/biosoftware/plink/v1.90beta4/plink --r2 inter-chr --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r2.ptt --allow-extra-chr
#/data/biosoftware/plink/v1.90beta4/plink --r inter-chr --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r.ptt --allow-extra-chr

#awk '$1!=$4  {print $0}' all.r.ptm.ld >all.r.ptm.inter.ld
#awk '$1!=$4  {print $0}' all.r2.ptm.ld >all.r2.ptm.inter.ld

#awk '$1!=$4  {print $0}' all.r.ptt.ld >all.r.ptt.inter.ld
#awk '$1!=$4  {print $0}' all.r2.ptt.ld >all.r2.ptt.inter.ld


####### Get Parental and Hybrid Frequencies
#vcftools --vcf all.recomblocks.ptm.recode.vcf --out all.recomblocks.ptm.recode --plink

#while read p 
#do
#  pairs=`echo "$p"| grep -v "SNP"|awk '{print $3 "\t" $6}'` 
#  /data/biosoftware/plink/v1.07/plink-1.07-x86_64/plink --ld $pairs --file all.recomblocks.ptm.recode --noweb \
#  |grep "AA\|TT\|AT\|TA\|A0\|T0\|0A\|0T" |awk '{print $2}' |grep -v "In\|LD\|ld\|info\|phase"| awk 'BEGIN { ORS = " " } { print }'
#  echo "$pairs"
#done <all.r.ptm.ld >all.pairs.freq

#/data/biosoftware/plink/v1.90beta4/plink --ld chrPTT_1:129482-129524 chrPTT_1:129524-131079 --vcf all.recomblocks.ptt.recode.vcf --recode structure --out test --allow-extra-chr 
#/data/biosoftware/plink/v1.90beta4/plink --ld chrPTT_1:129482-129524 0-1_contig_m86:100244-101644 --vcf all.recomblocks.ptt.recode.vcf --recode structure --out test2 --allow-extra-chr


###################################### Test Parallelization of Parental and Hybrid Frequencies f(x)
#mkdir splitfiles
#cd splitfiles
#cd splitfiles2
#cp ../all.r.ptm.inter.ld .
cp ../all.recomblocks.ptm.recode.vcf .

split -l 182400 -d all.r.ptm.inter.ld all.r.ptm.inter.ld-

for file in all.r.ptm.inter.ld-*
do

sbatch SL_recombination.v3.sh $file

done

echo "############## calculate the P-value on R for the top 0.01% ORs ####################"

cat all.r.ptm.inter.ld*txt |grep "chr.*chr\|chr.*m86\|m86.*chr"> all.odds
###at least 93 individuals represented
awk '{print $0 "\t" $3+$4+$5+$6}' all.odds|awk '$9 >= 93 {print $0}' >all.90indiv.odds 
###select top 0.05 % highest odds ratios
sort -V  -k7 all.90indiv.odds |awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}' |tail -7254 >all.01perc.odds

grep "chr\|contig" all.01perc.odds| awk '{print $1,$2,$0}'|awk  'sub(":","\t",$1)' | awk  'sub(":","\t",$3)' \
|awk  'sub("-","\t",$2)' | awk  'sub("-","\t",$5)'>all.tmp.odds

module load R/4.0.0
R
#install.packages("epitools")
library(epitools)
oddsratio<-read.table("all.tmp.odds")

p.value<-NULL

f <- function(x) {
 matrix<-matrix(as.numeric(as.character(x[c("V9","V10","V11","V12")])),ncol=2)
 oddsratio(matrix, method="fisher")$p.value[2,2]
}

p.value<-apply(oddsratio, 1, f)
oddsratio$p.value<-p.value
oddsratio$BH =p.adjust(oddsratio$p.value, method = "BH") ###FDR
oddsratio_sig<-subset(oddsratio, BH < 0.01,)
write.table(oddsratio_sig, file = "oddsratio_sigBH01.01perc.txt", quote = FALSE, sep="\t",row.names = FALSE,col.names=FALSE)
q()
###################### R LinkageDisequalibrium Heatmaps and top LD

#module load R
#Rscript --vanilla LinkageDisequilibrium.R all.r.ptm.matrix.ld all.r2.ptm.matrix.ld all.r2.ptm.ld all.r.ptm.inter.ld
#Rscript --vanilla LinkageDisequilibrium.R all.r.ptt.matrix.ld all.r2.ptt.matrix.ld all.r2.ptt.ld all.r.ptt.inter.ld

######## Identify R2 ranges that intersect with genic regions for both side of the pairings
sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n r_inter.1perc.trim.rpos.bed |cut -f1,4|uniq > r_inter_pairs

while read p;
do
  ##echo "$p"
  pairs=`echo "$p"|sed -e 's/\s\+/.*/g'`
  echo "$pairs"
  file=`echo "$p" |sed -e 's/\s\+/-/g'`
  grep ""$pairs"" r_inter.1perc.trim.rpos.bed |awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' \
  |bedtools merge -d 10000 -i stdin -o distinct,distinct,distinct,mean,mean,min,max -c 4,5,6,7,8,8,8 \
  |awk -v OFS='\t' '{gsub(",.*","",$5);print}' |awk -v OFS='\t' '{gsub(",.*","",$6);print}' \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  |bedtools intersect -a stdin -b all_genicregions.bed -wo |bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8,9,10 -c 14 -o distinct \
  |awk -v OFS="\t" '{print $4,$5,$6, $1,$2,$3, $7,$8, $11}' \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  |bedtools intersect -a stdin -b all_genicregions.bed -wo |bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8,9 -c 13 -o distinct >$file.r_inter.bed 
done <r_inter_pairs

cat chr*chr*r_inter.bed |awk -v OFS="\t" '{print $0,$7*$7}' |sort -k11,11 >r_inter.1perc.trim.rpos.genic.bed

awk -v OFS="\t" '{print $1,$4, $2,$3, $5,$6, $9,$10,$11}' r_inter.1perc.trim.rpos.genic.bed >r_inter.1perc.trim.rpos.genic.show-coords.txt




#########Identify Odds Ratio ranges that intersect with genic regions for both side of the pairings
sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n  oddsratio_sigBH01.01perc.txt > all.01perc.sort.odds
cat all.01perc.sort.odds|awk -v OFS="\t" '{print $1,$4}'|uniq >oddsratio_pairs

#####interatively collapse Odds Ratio ranges
while read p;
do
  ##echo "$p"
  pairs=`echo "$p"|sed -e 's/\s\+/.*/g'`
  echo "$pairs"
  file=`echo "$p" |sed -e 's/\s\+/-/g'`
  grep ""$pairs""  all.01perc.sort.odds|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' \
  |bedtools merge -d 10000 -i stdin -o distinct,min,max,mean,min,max -c 4,5,6,13,13,13 \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  |awk -v OFS="\t" '{print $4,$5,$6, $1,$2,$3, $7,$8,$9, $10}' \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  >$file.oddsratio.bed \

done <oddsratio_pairs

cat chr*chr*oddsratio.bed |awk -v OFS="\t" '{print $0,$7}' |awk 'BEGIN{OFS=FS="\t"}{ $12=sprintf("%.0f",$10) }1' |sort -k10,10n >odds_top.1perc.trim.bed

#### intersect with PTT genes ####
bedtools intersect -a odds_top.1perc.trim.bed -b ../all_genicregions.bed -wo |awk -v OFS="\t" '{print $0}'|bedtools groupby -i stdin -g 1,2,3,4,5,6 -c 16 -o distinct,count \
|awk -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8}' |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
|bedtools intersect -a stdin  -b ../all_genicregions.bed -wo |awk -v OFS="\t" '{print $0}'|bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o distinct,count > odds_top.1perc.trim.genes.bed

awk -v OFS="\t" '{print $1,$4, $2,$3, $5,$6, $7,$11}' odds_top.1perc.trim.genes.bed >odds_top.1perc.trim.show-coords.txt

#####get choords for PTM genome ######
awk -v OFS="\t" '{print $1,$3,$4,$2,$5,$6}' /home/yuzon/PTTxPTMtest/PTTxPTMcompare/circlize/PTTvsPTM.nuclearchrmitochondria.masked.show-coords.txt >PTTvsPTM.nuclearchrmitochondria.masked.show-coords.bed
bedtools sort -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai -i odds_top.1perc.trim.bed \
|bedtools intersect -a stdin -b PTTvsPTM.nuclearchrmitochondria.masked.show-coords.bed -wb -nonamecheck  \ 
###|awk  -v OFS="\t" '{gsub(/chrPTT_/, " ", $1)}1' |awk  -v OFS="\t" '{print $0,$15}'|awk  -v OFS="\t" '{gsub(/chrPTM_/, " ", $18)}1' |awk  -v OFS="\t" '$18 ==$1 {print $0}' comment out stringent chromosome match between PTT and PTM (e.g. chrPTT_3 == chrPTM_3 but chrPTT_3 !== chrPTM_2
|cut -f4-18 |bedtools sort -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai -i stdin \
|bedtools intersect -a stdin -b PTTvsPTM.nuclearchrmitochondria.masked.show-coords.bed -wb -nonamecheck |awk -v OFS="\t" '{print $12,$13,$14,$18,$19,$20}' \
|awk  -v OFS="\t" '$1 !=$4 {print $0}' >oddsratio_sigBH01.01perc.ptm.txt

sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n oddsratio_sigBH01.01perc.ptm.txt |awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t'  > all.01perc.sort.ptm.odds
cat all.01perc.sort.ptm.odds|awk -v OFS="\t" '{print $1,$4}'|uniq >oddsratio_pairs.ptm

#####interatively collapse Odds Ratio ranges
while read p;
do
  ##echo "$p"
  pairs=`echo "$p"|sed -e 's/\s\+/.*/g'`
  echo "$pairs"
  file=`echo "$p" |sed -e 's/\s\+/-/g'`
  grep ""$pairs""  all.01perc.sort.ptm.odds|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' \
  |bedtools merge -d 10000 -i stdin   -o distinct,min,max -c 4,5,6 \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  |awk -v OFS="\t" '{print $4,$5,$6, $1,$2,$3, $7,$8,$9, $10}' | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' \
  |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
  >$file.oddsratio.ptm.bed \

done <oddsratio_pairs.ptm

cat chr*chr*oddsratio.ptm.bed |awk -v OFS="\t" '{print $0,$7}' >odds_top.1perc.trim.ptm.bed


#### intersect with PTM	genes ####
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' odds_top.1perc.trim.ptm.bed |cut -f1-6 \
|bedtools sort -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai -i stdin| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' \
|bedtools intersect -a stdin -b ../all_genicregions.bed -wo |awk -v OFS="\t" '{print $0}'|bedtools groupby -i stdin -g 1,2,3,4,5,6 -c 10 -o distinct,count \
|awk -v OFS="\t" '{print $4,$5,$6,$1,$2,$3,$7,$8}' |bedtools sort -i stdin -faidx /home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta.fai \
|bedtools intersect -a stdin  -b ../all_genicregions.bed -wo |awk -v OFS="\t" '{print $0}' | sort -k1,1  -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n \
|bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o distinct,count	> odds_top.1perc.trim.genes.ptm.bed

awk -v OFS="\t" '{print $4,$1, $5,$6, $2,$3, $9,$7}' odds_top.1perc.trim.genes.ptm.bed >odds_top.1perc.trim.ptm.show-coords.txt

cat odds_top.1perc.trim.ptm.show-coords.txt odds_top.1perc.trim.show-coords.txt >odds_top.1perc.trim.ptmptt.show-coords.txt

#### get number of genes, effectors, and bsc
awk -v OFS="\t" '{print $1,$3,$4,$1 ":" $3 "-" $4 "_" $2 ":" $5 "-" $6}' odds_top.1perc.trim.ptmptt.show-coords.txt >ptt_first
awk -v OFS="\t" '{print $2,$5,$6,$1 ":" $3 "-" $4 "_" $2 ":" $5 "-" $6}' odds_top.1perc.trim.ptmptt.show-coords.txt >ptt_second 
cat ptt_first ptt_second | bedtools intersect -a stdin -b ../all_genicregions.bed -wo | sort -k4,4 |bedtools groupby -g 4 -c 8 -o freqasc |sort >odds_top.1perc.trim.ptmptt.counts
### removed chrPTT_1:2520786-2522920_chrPTT_7:1932789-1932834 because no genes in the PTT genome therefore no DMI

