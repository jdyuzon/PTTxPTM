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
#SBATCH --time=10:00:00
#  maximum requested memory
#SBATCH --mem=15G  #70G
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

date
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
#bwa index ${main_dir%/Populations/}/references/hybrid_references/merged_genome_PTTxPTM.fasta

Raw_reads=${main_dir/Populations/}/raw_reads/Iranianpop/ #PTTpop/
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

gridss=${main_dir}gridss
manta=${main_dir}manta
parsnp_qry=${main_dir}SNPcalling/parsnp_qry/
parsnp_out=${main_dir}SNPcalling/parsnp_out/
nucmer=${main_dir}SNPcalling/nucmer/

mkdir $refPTT
mkdir $refPTM
mkdir $merged_genome

mkdir $Raw_reads
mkdir $Raw_fastqc
mkdir $trimmed

mkdir $PEoverlap
mkdir $Pear_fastqc
mkdir $assembled
mkdir $unassembled
mkdir $unassembled_fastqc
mkdir $assembled_fastqc

mkdir $markdup
mkdir $markdup

mkdir $bin_genome

mkdir $readgroup_files

mkdir $bwa_out

mkdir $stampy_realign
mkdir $flag_out

mkdir $Sorted
mkdir $Depth_dir

mkdir $snps

mkdir $gridss
mkdir $manta
mkdir $parsnp_qry
mkdir $parsnp_out

mkdir $nucmer

#####################################################
echo "Trim adaptors, min quality =30, min length =40"
#####################################################
cd $Raw_reads
echo $Raw_reads
for each in ${entry}_R1.fastq.gz
do
echo ${each}
java -jar ${trimmo}trimmomatic-0.38.jar PE -threads 5 \
      $each ${entry}_R2.fastq.gz \
      ${trimmed}${entry}_P_1.fastq ${trimmed}${entry}_U_1.fastq \
      ${trimmed}${entry}_P_2.fastq ${trimmed}${entry}_U_2.fastq \
      ILLUMINACLIP:${trimmo}adapters/TruSeq3-PE-2.fa:2:30:10 \
      SLIDINGWINDOW:3:28 MINLEN:40
done

##################FastQC##############################

cd ${trimmed}
${fastqc} ${entry}_P_1.fastq
${fastqc} ${entry}_P_2.fastq
${fastqc} ${entry}_U_1.fastq
${fastqc} ${entry}_U_2.fastq

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


##########################################################
echo "BWA: Alignment"#####################################
##########################################################
cd ${markdup}

##${bwa_soft}bwa index ${refPTT}
##${bwa_soft}bwa index ${refPTM}

#####map to PTT
for each in ${entry}.markdup_1.fastq
do
isolate=${each%.markdup_1.fastq}
echo $isolate
${bwa_soft}bwa mem -M -t 1  ${refPTT} ${each} ${isolate}.markdup_2.fastq \
> ${bwa_out}${isolate}.PTT.PE.sam

${bwa_soft}bwa mem -M -t 1 ${refPTT} \
${isolate}.markdup.fastq > ${bwa_out}${isolate}.PTT.SE.sam
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

rm ${entry}.ptg*
${stampy_soft}stampy.py -G ${entry}.ptg ${refPTT}
${stampy_soft}stampy.py -g ${entry}.ptg -H ${entry}.ptg

for each in ${entry}.PTT.SE.bam
do
echo ${each}
${stampy_soft}stampy.py -g ${entry}.ptg -h ${entry}.ptg -t4 --bamkeepgoodreads -M ${each} -o ${bwa_out}${each%.bam}.stampy.sam --overwrite
${stampy_soft}stampy.py -g ${entry}.ptg -h ${entry}.ptg -t4 --bamkeepgoodreads -M ${each%.SE.bam}.PE.bam -o ${bwa_out}${each%.SE.bam}.PE.stampy.sam --overwrite
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
## stampy depth stats
samtools depth -a ${Sorted}${each%.stampy.sam}.sorted.bam > ${Depth_dir}${each}.stats
done

##########################################################
echo "Samtools: Merge"####################################
##########################################################
cd ${Sorted}

samtools merge -f ${entry}.PTMPTT.merged.bam \
$entry.PTT.PE.sorted.bam $entry.PTT.SE.sorted.bam \

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
       RGID=${entry} \
       RGLB=${entry} \
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
   --annotation MappingQuality \
   -O ${snps}${each%.read.bam}.vcf.gz

done


##########################################################
echo "BCFtools: set low AD to missing "###################
##########################################################
cd $snps
###set indels to missing

bcftools filter -e 'FORMAT/AD[0] < 5 & FORMAT/GT = "0" | FORMAT/AD[1] < 5 & FORMAT/GT = "1" | FORMAT/AD[2] < 5 & FORMAT/GT = "2"'  --set-GTs . ${entry}.PTMPTT.vcf.gz \
| bcftools filter -e ' ((FORMAT/AD[0])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "0" |  ((FORMAT/AD[1])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "1"|  ((FORMAT/AD[2])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "2"'  --set-GTs . \
| bcftools filter -e ' %QUAL<20 ' --set-GTs . | bcftools filter -e 'TYPE="indel" | TYPE="mnp" | TYPE="other" '  --set-GTs . | sed 's/	\.:/	\.\/.:/g'> ${entry}.PTMPTT.flt.vcf

##########################################################
echo "Bcftools consensus: VCF to FASTA "##################
##########################################################
bgzip -f -c ${entry}.PTMPTT.flt.vcf > ${entry}.PTMPTT.flt.vcf.gz
tabix -f -p vcf ${entry}.PTMPTT.flt.vcf.gz
####for variant sites remove ",<NON_REF>" so only ALT nucleotide, then for sites with no calls replace "<NON_REF>" with "N". ####
####Double check that REF/ALT only has SNPs no indels. This will remove all indel sites #####
cat ${entry}.PTMPTT.flt.vcf |sed 's/,<NON_REF>//g' | sed 's/<NON_REF>/N/g' |sed 's/\*/N/g' \
|awk 'length($4)==1 || $1 ~/#/' | awk 'length($5)==1|| $1 ~/#/' \
|bgzip -c -f > ${entry}.PTMPTT.flt2.vcf.gz
tabix -f -p vcf ${entry}.PTMPTT.flt2.vcf.gz


#### Create consensus fasta and mask sites ####
/data/biosoftware/bcftools/bcftools-1.11/bcftools consensus -c $entry.chain -f $refPTT ${entry}.PTMPTT.flt2.vcf.gz -H 1 \
|sed 's/\*/N/g'  > ${entry}.PTMPTT.fasta



############################################################
#### Get all genes from Fasta ########################### ##
############################################################

##remove old index file
#rm ${entry}.PTMPTT.fasta.fai

#### make sure that you get the right bed ranges because vcf to fasta can change the coordinates
bedtools getfasta -tab -fi ${entry}.PTMPTT.fasta -bed ~/PTTxPTMtest/PTTxPTMcompare/SegregationDistortion_files/paml/paml_pttptm_2/allgenes.bed \
|sed 's/N/?/g'| awk -v variable="$entry" '{print $1 "_" variable "\t" $2} '  >${entry}.allgenes.fasta

###replace N with ? (undetermined nucleotide or amino acid) for PAML 

############################################################
#### Bedtools: flatten vcf to bed ranges ######
############################################################
awk '{print $1 "\t" $2 "\t" $2 "\t" $10 "\t" $10}' ${snps}${entry}.PTMPTT.flt.vcf \
|grep -v "#"|sed 's/:.*//g'|awk -v OFS='\t' '{ if ($4 == "0") {$4 = "+"}; print }'\
| awk -v OFS='\t'  '{ if ($4 == "1") {$4 = "-"}; print }' |awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $3 "\t" $3 "\t"$4}'>${entry}.tmp.bed
bedtools merge -i ${entry}.tmp.bed -s -c 6 -o distinct | awk 'BEGIN { OFS = "\t" } { $5 = $2 - $3 } 1' \
|awk -v variable="$entry" '{print $0 "\t" variable} ' >${entry}.PTMPTT.flt.bed


############################################################
#### Nucmer: get PTT coordinates for recombination blocks ##
############################################################
cd ${nucmer}

nucmer -mum -mincluster 100 -minmatch 1000 --prefix=${nucmer}${entry} ~/references/hybrid_references/PTT_0-1_assembly.v14.fa ${snps}${entry}.PTMPTT.fasta
delta-filter -r -q ${nucmer}${entry}.delta >${nucmer}${entry}.filter
show-snps -Clr ${nucmer}${entry}.filter > ${nucmer}${entry}.snps
show-coords ${nucmer}${entry}.delta > ${nucmer}${entry}.show-coords
grep -v "contig.*PTM" ${nucmer}${entry}.show-coords \
|awk '{print $12 "\t" $13 "\t" $1 "\t" $2 "\t" $4 "\t" $5}'>${nucmer}${entry}.show-coords.txt


date

