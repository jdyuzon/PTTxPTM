#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=PTTxPTMhybrid
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=12:00:00
#  maximum requested memory
#SBATCH --mem=15G #70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTMhybrid.%J.err
#SBATCH --output=PTTxPTMhybrid.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

###########################################################
echo "Input Isolate Name" $entry###########################
###########################################################
### example: sbatch SL_PTTxPTMhybrid.parallel.sh 100 YOUR_HOME_DIR 
entry=$1 
main_dir=$2

echo $entry
###########################################################
echo "Software"
###########################################################
##module load java/x64/8u121

trimmo=/data/biosoftware/Trimmomatic/Trimmomatic-0.38/
pear=/data/biosoftware/pear/pear-0.9.10-bin-64/pear-0.9.10-bin-64 
export BBMAPDIR=/data/biosoftware/bbmap/bbmap
bwa_soft=/data/biosoftware/bwa/bwa-0.7.15/
stampy_soft=/data/biosoftware/stampy/stampy/
samtools_soft=/data/biosoftware/samtools/samtools-1.9/
fastqc=/data/biosoftware/FastQC/FastQC/fastqc


##########################################################
echo "Directories"
##########################################################
refPTT=~/references/hybrid_references/PTT_0-1_assembly.v14.fa
refPTM=~/references/hybrid_references/PTM_FGOB10Ptm-1_assembly.v7.fasta
merged_genome=~/references/hybrid_references/merged_genome_PTTxPTM.fasta
####index fasta
#bwa index ${main_dir%/PTTxPTMtest/}/references/hybrid_references/merged_genome_PTTxPTM.fasta

Raw_reads=${main_dir/PTTxPTMtest/}/raw_reads/
Raw_fastqc=${main_dir}fastQC/raw_fastqc
trimmed=${main_dir}trimmed/

PEoverlap=${main_dir}trimmed/pear/
assembled=${main_dir}trimmed/pear/assembled/
unassembled=${main_dir}trimmed/pear/unassembled/
unassembled_fastqc=${main_dir}fastQC/fastqc_unassembled
assembled_fastqc=${main_dir}fastQC/fastqc_assembled

markdup=${main_dir}trimmed/prinseq/
markdup_fastQC=${main_dir}fastQC/markdup

bin_genome=${main_dir}bbsplit/
bwa_out=${main_dir}aligned/

stampy_realign=${main_dir}aligned/stampy/
flag_out=${main_dir}stampy/flag/
readgroup_files=${main_dir}aligned/sorted/readgroups/

Sorted=${main_dir}aligned/sorted/
Depth_dir=${main_dir}aligned/sorted/stats/depth/
dups_dir=${main_dir}aligned/sorted/stats/dups/

snps=${main_dir}SNPcalling/

parsnp_qry=${main_dir}SNPcalling/parsnp_qry/
parsnp_out=${main_dir}SNPcalling/parsnp_out/
nucmer=${main_dir}SNPcalling/nucmer/
nucmer_PTM=${main_dir}SNPcalling/nucmer_PTM/

#mkdir $refPTT
#mkdir $refPTM
#mkdir $merged_genome

#mkdir $Raw_reads
#mkdir $Raw_fastqc
#mkdir $trimmed

#mkdir $PEoverlap
#mkdir $Pear_fastqc
#mkdir $assembled
#mkdir $unassembled
#mkdir $unassembled_fastqc
#mkdir $assembled_fastqc

#mkdir $markdup
#mkdir $markdup

#mkdir $bin_genome

#mkdir $readgroup_files

#mkdir $bwa_out

#mkdir $stampy_realign
#mkdir $flag_out

#mkdir $Sorted
#mkdir $Depth_dir

#mkdir $snps

#mkdir $gridss
#mkdir $manta
#mkdir $parsnp_qry
#mkdir $parsnp_out

mkdir $nucmer
mkdir $nucmer_PTM
#####################################################
echo "Trim adaptors, min quality =30, min length =40"
#####################################################
cd ${Raw_reads}
for each in ${entry}_1.fq.gz
do
echo ${each}
java -jar ${trimmo}trimmomatic-0.38.jar PE -threads 5 \
      $each ${each%1.fq.gz}2.fq.gz \
      ${trimmed}${each%1.fq.gz}P_1.fastq ${trimmed}${each%1.fq.gz}U_1.fastq \
      ${trimmed}${each%1.fq.gz}P_2.fastq ${trimmed}${each%1.fq.gz}U_2.fastq \
      ILLUMINACLIP:${trimmo}adapters/TruSeq3-PE-2.fa:2:30:10 \
      SLIDINGWINDOW:3:28 MINLEN:40
done

##################FastQC##############################

#cd ${trimmed}
#${fastqc} ${entry}_P_1.fastq
#${fastqc} ${entry}_P_2.fastq
#${fastqc} ${entry}_U_1.fastq
#${fastqc} ${entry}_U_2.fastq

##########################################################################################################
echo "Merge Singleton Files: Trimmomatic Unpaired (For and Rev) " ########################################
##########################################################################################################
for each in ${entry}_U_1.fastq
do
cat ${each} \
${each%_U_1.fastq}_U_2.fastq \
> ${assembled}${each%_U_1.fastq}.singles.fastq

done

###########################################################
echo "Prinseq: remove PCR duplicates"######################
###########################################################
###prinseq needs this version of perl but this version interferes with other software

module load perl/5.26.1

cd ${assembled}
################## singles #################################
for each in ${entry}.singles.fastq
do
echo ${each}
perl /data/biosoftware/prinseq/bin/prinseq-lite.pl -fastq ${each} \
-log -verbose -derep 14 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 30 -trim_qual_step 3 \
-min_len 40 -out_format 3 -out_good ${markdup}${each%.singles.fastq}.markdup
done

##paired
cd ${trimmed}

for each in ${entry}_P_1.fastq 
do 
echo ${each}
echo ${each%_P_1}_P_2.fastq 
perl /data/biosoftware/prinseq/bin/prinseq-lite.pl \
-fastq ${each} \
-fastq2 ${each%_P_1.fastq}_P_2.fastq \
-log -verbose -derep 14 -out_format 3 \
-out_good ${markdup}${each%_P_1.fastq}.markdup

done

###prinseq needs this version of perl but this version interferes with other software
module unload perl/5.26.1 

for each in ${entry}_P_1.fastq
do
fastqc ${markdup}${each%_P_1.fastq}.markdup_1.fastq --outdir ${markdup_fastqc}
fastqc ${markdup}${each%_P_1.fastq}.markdup_2.fastq --outdir ${markdup_fastqc}
fastqc ${markdup}${each%_P_1.fastq}.markdup.fastq --outdir ${markdup_fastqc}
done


#################################################################
echo "BBSplit: Bin reads to best matching Parent Genome"#########
#################################################################
cd ${markdup}

for each in ${entry}.markdup_1.fastq
do
#####paired (ambiguous2=all map to both genomes)
$BBMAPDIR/bbsplit.sh minratio=0.52 ambiguous=best ambiguous2=all \
in1=${each} in2=${each%.markdup_1.fastq}.markdup_2.fastq \
ref=$refPTT,$refPTM basename=${bin_genome}${each%.markdup_1.fastq}_P_%.fq
#####singles
$BBMAPDIR/bbsplit.sh minratio=0.52 ambiguous=best ambiguous2=all \
in1=${each%.markdup_1.fastq}.markdup.fastq \
ref=$refPTT,$refPTM basename=${bin_genome}${each%.markdup_1.fastq}_%.fq
done

##########################################################
echo "BWA: Alignment"#####################################
##########################################################
cd ${bin_genome}

${bwa_soft}bwa index ${refPTT}
${bwa_soft}bwa index ${refPTM}

#####map to PTT
for each in ${entry}_P_PTT_0-1_assembly.v14.fq
do
isolate=${each%_P_PTT_0-1_assembly.v14.fq}
echo $isolate
${bwa_soft}bwa mem -M -t 4 -p ${refPTT} ${each} \
> ${bwa_out}${each%_P_PTT_0-1_assembly.v14.fq}.PTT.PE.sam

${bwa_soft}bwa mem -M -t 4 ${refPTT} \
${isolate}_PTT_0-1_assembly.v14.fq > ${bwa_out}${isolate}.PTT.SE.sam
done


#####map to PTM
for each in ${entry}_P_PTM_FGOB10Ptm-1_assembly.v7.fq
do
isolate=${each%_P_PTM_FGOB10Ptm-1_assembly.v7.fq}
echo $isolate
${bwa_soft}bwa mem -M -t 4 -p ${refPTM} ${each} \
> ${bwa_out}${each%_P_PTM_FGOB10Ptm-1_assembly.v7.fq}.PTM.PE.sam

${bwa_soft}bwa mem -M -t 4 ${refPTM} \
${isolate}_PTM_FGOB10Ptm-1_assembly.v7.fq > ${bwa_out}${isolate}.PTM.SE.sam

done



##########################################################
echo "Stampy: Indel Realigner/Assess Insert Size"#########
##########################################################
module load python/2.7.13

cd ${bwa_out}
####make BAM
for each in ${entry}*SE.sam
do
echo ${each}
${samtools_soft}samtools view -S -b ${each} > ${each%.sam}.bam  
${samtools_soft}samtools view -S -b ${each%.SE.sam}.PE.sam > ${each%.SE.sam}.PE.bam
done

####stampy
rm ptg*
${stampy_soft}stampy.py -G ptg ${refPTT}
${stampy_soft}stampy.py -g ptg -H ptg

for each in ${entry}.PTT.SE.bam
do
echo ${each}
${stampy_soft}stampy.py -g ptg -h ptg -t4 --bamkeepgoodreads -M ${each} -o ${bwa_out}${each%.bam}.stampy.sam
${stampy_soft}stampy.py -g ptg -h ptg -t4 --bamkeepgoodreads -M ${each%.SE.bam}.PE.bam -o ${bwa_out}${each%.SE.bam}.PE.stampy.sam
done

rm ptg*
${stampy_soft}stampy.py -G ptg ${refPTM}
${stampy_soft}stampy.py -g ptg -H ptg

for each in ${entry}.PTM.SE.bam
do
echo ${each}
${stampy_soft}stampy.py -g ptg -h ptg -t4 --bamkeepgoodreads -M ${each} -o ${bwa_out}${each%.bam}.stampy.sam
${stampy_soft}stampy.py -g ptg -h ptg -t4 --bamkeepgoodreads -M ${each%.SE.bam}.PE.bam -o ${bwa_out}${each%.SE.bam}.PE.stampy.sam
done

cd ${bwa_out}
for each in ${entry}*stampy.sam
do
echo ${each}
${samtools_soft}samtools flagstat ${each} > ${flag_out}${each}.stats
done

module unload python/2.7.13

##########################################################
echo "Samtools: Sort"#####################################
##########################################################
cd ${bwa_out}
for each in ${entry}*stampy.sam
do
echo ${each}
samtools sort -o ${Sorted}${each%.stampy.sam}.sorted.bam -@ 4 ${each}
### stampy depth stats
samtools depth -a ${Sorted}${each%.stampy.sam}.sorted.bam > ${Depth_dir}${each}.stats
done

##########################################################
echo "Samtools: Merge"####################################
##########################################################
cd ${Sorted}

for each in ${entry}.PTM.PE.sorted.bam
do
echo ${each}
samtools merge ${each%.PTM.PE.sorted.bam}.PTMPTT.merged.bam \
${each} ${each%.PE.sorted.bam}.SE.sorted.bam \
${each%.PTM.PE.sorted.bam}.PTT.PE.sorted.bam ${each%.PTM.PE.sorted.bam}.PTT.SE.sorted.bam
done

##########################################################
echo "Picard: ReadGroup"##################################
##########################################################
for each in ${entry}.PTMPTT.merged.bam
do
echo ${each}
var=`samtools view ${each} | head -n 1 | cut -f1 | cut -d ":" -f 3,4`
java -jar /data/biosoftware/Picard/Picard/picard.jar AddOrReplaceReadGroups \
       I=${each} \
       O=${each%.merged.bam}.read.bam \
       RGID="$var" \
       RGLB=lb1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=${each%.merged.bam}
samtools index ${each%.merged.bam}.read.bam
done

##########################################################
echo "GATK: Haplotype Caller"#############################
##########################################################
cd ${Sorted}
for each in ${entry}.PTMPTT.read.bam
do

samtools index $each
/data/biosoftware/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx4g" HaplotypeCaller \
   -R ${merged_genome} \
   -I ${Sorted}${each} \
   -ERC BP_RESOLUTION \
   -ploidy 1 \
   -O ${snps}${each%.read.bam}.vcf.gz

done


##########################################################
echo "BCFtools: set low AD to missing "###################
##########################################################
cd $snps

bcftools filter -e 'FORMAT/AD[*:0] < 3 | FORMAT/AD[*:1] < 3 | FORMAT/AD[*:2] < 3' --set-GTs .  ${entry}.PTMPTT.vcf.gz > ${entry}.PTMPTT.flt.vcf
##########################################################
echo "Bcftools consensus: VCF to FASTA "##################
##########################################################

grep -v "#" ${entry}.PTMPTT.flt.vcf |awk '{print $1 "\t" $2 "\t" $2 "\t" $10}'|sed 's/:.*//g'|awk '$4 ~ /\./||$4 ~ /1/ || $4 ~ /2/' > ${entry}.mask.bed
bedtools maskfasta -fi $merged_genome -bed ${entry}.mask.bed -fo ${entry}.PTMPTT.fasta

############################################################
#### Split isolate.fasta files to PTT and PTM ##############
############################################################
bedtools getfasta  -fi ${entry}.PTMPTT.fasta -bed ~/references/hybrid_references/PTM_tracks.bed|sed 's/:.*//g' > ${parsnp_qry}${entry}.PTM.fasta
bedtools getfasta  -fi ${entry}.PTMPTT.fasta -bed ~/references/hybrid_references/PTT_tracks.bed|sed 's/:.*//g' > ${parsnp_qry}${entry}.PTT.fasta



############################################################
#### Bedtools: flatten vcf to bed ranges ######
############################################################
awk '{print $1 "\t" $2 "\t" $2 "\t" $10 "\t" $10}' ${snps}${entry}.PTMPTT.flt.vcf \
|grep -v "#"|sed 's/:.*//g'|awk -v OFS='\t' '{ if ($4 == "0") {$4 = "+"}; print }'\
| awk -v OFS='\t'  '{ if ($4 == ".") {$4 = "-"}; print }' |awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $3 "\t" $3 "\t"$4}'>${entry}.tmp.bed
bedtools merge -i ${entry}.tmp.bed -s -c 6 -o distinct | awk 'BEGIN { OFS = "\t" } { $5 = $2 - $3 } 1' \
|awk -v variable="$entry" '{print $0 "\t" variable} ' >${entry}.PTMPTT.flt.bed


############################################################
#### Nucmer: get PTT coordinates for recombination blocks ##
############################################################
nucmer -mum -mincluster 100 -minmatch 1000 --prefix=${nucmer}${entry} ~/references/hybrid_references/PTT_0-1_assembly.v14.fa ${entry}.PTMPTT.fasta
delta-filter -r -q ${nucmer}${entry}.delta >${nucmer}${entry}.filter
show-snps -Clr ${nucmer}${entry}.filter > ${nucmer}${entry}.snps
show-coords ${nucmer}${entry}.delta > ${nucmer}${entry}.show-coords
grep -v "contig.*PTM" ${nucmer}${entry}.show-coords \
|awk '{print $12 "\t" $13 "\t" $1 "\t" $2 "\t" $4 "\t" $5}'>${nucmer}${entry}.show-coords.txt

nucmer -mum -mincluster 65 -minmatch 40 --prefix=${nucmer}${entry}.alt  $refPTM ${entry}.PTMPTT.fasta
delta-filter -r -q ${nucmer}${entry}.alt.delta >${nucmer}${entry}.alt.filter
show-snps -Clr ${nucmer}${entry}.alt.filter > ${nucmer}${entry}.alt.snps

cd ${nucmer}

grep "chrPTM\|mito" ${entry}.show-coords|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.PTM.bed
grep -v "PTM\|mito" ${entry}.show-coords|grep "chrPTT\|contig"|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' |bedtools sort -i stdin -faidx ${refPTT}.fai>${entry}.PTT.bed

bedtools intersect -a ${entry}.PTT.bed -b ${entry}.PTM.bed -nonamecheck|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "common"}'> ${entry}.common.bed
bedtools subtract -a ${entry}.PTM.bed -b ${entry}.PTT.bed -nonamecheck|cut -f1,2,3|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "PTM"}'> ${entry}.PTMuniq.bed
bedtools subtract -a ${entry}.PTT.bed -b ${entry}.PTM.bed -nonamecheck|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "PTT"}'> ${entry}.PTTuniq.bed

cat ${entry}.common.bed ${entry}.PTMuniq.bed ${entry}.PTTuniq.bed|bedtools sort -i stdin -faidx ${refPTT}.fai >${entry}.common_uniq.bed

echo "###merge bed ranges: if common is between PTT alleles, change common to PTT; if common is between PTM alleles, change common to PTM###"
echo '.	-1	-1	.' | cat -  ${entry}.common_uniq.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_uniq.bed >${entry}.down.tmp
paste ${entry}.common_uniq.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp

awk -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTT" && $12=="PTT") $4="PTT" ; print}' ${entry}.updown.tmp \
| awk -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTM" && $12=="PTM") $4="PTM" ; print}' >${entry}.tmp

awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk -v OFS='\t' '$4 == "PTM"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp
cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "##second round of merging: if downstream allele same as current allele, change current STOP to downstream STOP; if upstream allele same as current allele, change current START to downstream START###"
echo '.	-1	 -1	.' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp

awk -v OFS='\t' '{if($1==$9 && $4==$12) $3=$11; print}' ${entry}.updown.tmp \
|awk -v OFS='\t' '{if($1==$5 && $4==$8) $2=$6; print}' \
|awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'|awk -v OFS='\t' '{print $0, "PTM"}' |bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'|awk -v OFS='\t' '{print $0, "PTT"}' |bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'|awk -v OFS='\t' '{print $0, "common"}' |bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp

cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp >${entry}.tmp
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' ${entry}.tmp \
|bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "###third round of merging: if common is between PTT alleles, change common to PTT; if common is between PTM alleles, change common to PTM###"
echo '.	-1	-1	.' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp
awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTT" && $12=="PTT") $4="PTT" ; print}' ${entry}.updown.tmp \
| awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTM" && $12=="PTM") $4="PTM" ; print}' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp
cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "##### Number of SNPs supporting block: chr start stop PTT/PTM isolate SNPs #####"
awk -v OFS='\t' '{print $14,$1,$1,$15}' ${nucmer}${entry}.snps |grep "chrPTM\|mito" |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.PTM.snps
awk -v OFS='\t' '{print $15,$4,$4,$15}' ${entry}.alt.snps |grep "chrPTT\|0-1_contig" |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.PTT.snps
######grep  "chrPTT\|0-1_contig_" $snps${entry}.PTMPTT.flt.vcf |grep -v "#" |awk '{print $1 "\t" $2 "\t" $2 "\t" $10}'|sed 's/:.*//g'| awk '$4 ~ /0/' |bedtools sort -i stdin  -faidx ${refPTT}.fai >${entry}.PTT.snps
bedtools subtract -a ${entry}.PTM.snps -b ${entry}.PTT.snps -nonamecheck |bedtools intersect -a $entry.ptm.tmp -b stdin -wa -nonamecheck |bedtools groupby -i stdin -c 1 -o count -nonamecheck |awk '$4 > 5' | awk -v variable="$entry" '{print $1 "\t" $2 "\t" $3 "\t" "PTM" "\t" variable "\t" $4}'|bedtools sort -i stdin -faidx ${refPTT}.fai  >${entry}.PTM.support
bedtools subtract -a ${entry}.PTT.snps -b ${entry}.PTM.snps -nonamecheck |bedtools intersect -a $entry.ptt.tmp -b stdin -wa -nonamecheck |bedtools groupby -i stdin -c 1 -o count -nonamecheck |awk '$4 > 5' | awk -v variable="$entry" '{print $1 "\t" $2 "\t" $3 "\t" "PTT" "\t" variable "\t" $4}'|bedtools sort -i stdin -faidx ${refPTT}.fai  >${entry}.PTT.support
awk -v variable="$entry" '{print $0 "\t"  variable}' $entry.common.tmp > tmp && mv tmp $entry.common.tmp

cat ${entry}.PTM.support ${entry}.PTT.support $entry.common.tmp|cut -f1-5|bedtools sort -i stdin -faidx ${refPTT}.fai  > ${entry}.common_unique.blocks.bed

echo "###fourth round of merging: if common is between PTT alleles, change common to PTT; if common is between PTM alleles, change common to PTM###"
echo '. -1	-1	.' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp
awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTT" && $12=="PTT") $4="PTT" ; print}' ${entry}.updown.tmp \
| awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTM" && $12=="PTM") $4="PTM" ; print}' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp
cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "##fifth round of merging: if downstream allele same as current allele, change current STOP to downstream STOP; if upstream allele same as current allele, change current START to downstream START######"
echo '. -1	 -1     .' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp

awk -v OFS='\t' '{if($1==$9 && $4==$12) $3=$11; print}' ${entry}.updown.tmp \
|awk -v OFS='\t' '{if($1==$5 && $4==$8) $2=$6; print}' \
|awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'|awk -v OFS='\t' '{print $0, "PTM"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'|awk -v OFS='\t' '{print $0, "PTT"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'|awk -v OFS='\t' '{print $0, "common"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp

cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp >${entry}.tmp
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' ${entry}.tmp \
|bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "###sixth round of merging: if common is between PTT alleles, change common to PTT; if common is between PTM alleles, change common to PTM###"
echo '. -1	-1	.' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp
awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTT" && $12=="PTT") $4="PTT" ; print}' ${entry}.updown.tmp \
| awk  -v OFS='\t' '{if($1==$5 && $1==$9 && $4=="common" && $8=="PTM" && $12=="PTM") $4="PTM" ; print}' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'| bedtools sort -i stdin -faidx ${refPTT}.fai \
| bedtools merge -d 1 -i stdin |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp
cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "##seventh round of merging: if downstream allele same as current allele, change current STOP to downstream STOP; if upstream allele same as current allele, change current START to downstream START####"
echo '. -1	 -1     .' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp > ${entry}.updown.tmp

awk -v OFS='\t' '{if($1==$9 && $4==$12) $3=$11; print}' ${entry}.updown.tmp \
|awk -v OFS='\t' '{if($1==$5 && $4==$8) $2=$6; print}' \
|awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' >${entry}.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'|awk -v OFS='\t' '{print $0, "PTM"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'|awk -v OFS='\t' '{print $0, "PTT"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'|awk -v OFS='\t' '{print $0, "common"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp

cat $entry.ptm.tmp $entry.ptt.tmp $entry.common.tmp >${entry}.tmp
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' ${entry}.tmp \
|bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.common_unique.blocks.bed

echo "###################### Recombination Breakpoint #############################"
echo '.	-1	-1	.' | cat -     ${entry}.common_unique.blocks.bed >${entry}.up.tmp
tail -n +2 ${entry}.common_unique.blocks.bed >${entry}.down.tmp
paste ${entry}.common_unique.blocks.bed ${entry}.up.tmp ${entry}.down.tmp|awk -v OFS="\t" '$1=$1' > ${entry}.updown.tmp
awk -v OFS='\t' '{if($1==$9 && $1==$5 && $4=="common" && $8=="PTM" && $12=="PTT"){print $1,$2,$3,$4} \
	     else if($1==$9 && $1==$5 && $4=="common" && $8=="PTT" && $12=="PTM"){print $1,$2,$3,$4} \
	     else if($1==$5 && $4=="PTM" && $8=="PTT"){print $1, $7, $2, "break"} \
	     else if($1==$5 && $4=="PTT" && $8=="PTM"){print $1, $7, $2, "break"}}' ${entry}.updown.tmp \
|awk -v variable="$entry" '{print $0 "\t" variable }'>$entry.break
	
	
