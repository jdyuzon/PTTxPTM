#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=PTTxPTTcross
#  how many cpus are requested
#SBATCH --ntasks=4
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=15G #70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTTcross.%J.err
#SBATCH --output=PTTxPTTcross.%J.out
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
refPTT=/home/yuzon/references/hybrid_references/PTT_0-1_assembly.v14.fa
refPTM=/home/yuzon/references/Wyatt_Ellwood_assemblies/15A.fasta
merged_genome=/home/yuzon/references/PTTxPTT_genome/merged_genome_PTTxPTT.fasta
####index fasta
#bwa index ${main_dir%/15Ax0-1_ptt/}/references/hybrid_references/merged_genome_PTTxPTT.fasta

Raw_reads=${main_dir/15Ax0-1_ptt/}/raw_reads/
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
nucmer_PTM=${main_dir}SNPcalling/nucmer_PTM/

#mkdir $refPTT
#mkdir $refPTM
#mkdir $merged_genome

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
mkdir $nucmer_PTM
#####################################################
echo "Trim adaptors, min quality =30, min length =40"
#####################################################
cd ${Raw_reads}

/data/biosoftware/bbmap/bbmap/bbduk.sh -Xmx1g \
in=$entry.fastq out=${trimmed}$entry.trim.fastq ref=/home/yuzon/BioSoft/IonTorrent_adapters.fa \
 ktrim=r k=23 mink=11 hdist=1 tpe tbo minlen=40 trimq=10 qtrim=rl

##################FastQC##############################

cd ${trimmed}
${fastqc} ${entry}.trim.fastq


#################################################################
echo "BBSplit: Bin reads to best matching Parent Genome"#########
#################################################################
cd ${trimmed}

for each in ${entry}.trim.fastq
do
#####singles
$BBMAPDIR/bbsplit.sh minratio=0.52 ambiguous=best ambiguous2=all \
in1=$each \
ref=$refPTT,$refPTM basename=${bin_genome}${each%.trim.fastq}_%.fq
done

##########################################################
echo "BWA: Alignment"#####################################
##########################################################
cd ${bin_genome}

###${bwa_soft}bwa index ${refPTT}
###${bwa_soft}bwa index ${refPTM}

#####map to PTT
${bwa_soft}bwa mem -M -t 4 ${refPTT} \
${entry}_PTT_0-1_assembly.v14.fq > ${bwa_out}${entry}.PTT.SE.sam

#####map to PTM
${bwa_soft}bwa mem -M -t 4 ${refPTM} \
${entry}_15A.fq > ${bwa_out}${entry}.PTM.SE.sam



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
done

####stampy
rm $entry.ptg*
${stampy_soft}stampy.py -G $entry.ptg ${refPTT}
${stampy_soft}stampy.py -g $entry.ptg -H $entry.ptg

for each in ${entry}.PTT.SE.bam
do
echo ${each}
${stampy_soft}stampy.py --overwrite -g $entry.ptg -h $entry.ptg -t4 --bamkeepgoodreads -M ${each} -o ${bwa_out}${each%.bam}.stampy.sam
done

rm $entry.ptg*
${stampy_soft}stampy.py -G $entry.ptg ${refPTM}
${stampy_soft}stampy.py -g $entry.ptg -H $entry.ptg

for each in ${entry}.PTM.SE.bam
do
echo ${each}
${stampy_soft}stampy.py --overwrite -g $entry.ptg -h $entry.ptg -t4 --bamkeepgoodreads -M ${each} -o ${bwa_out}${each%.bam}.stampy.sam
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

for each in ${entry}.PTM.SE.sorted.bam
do
echo ${each}
samtools merge ${entry}.PTMPTT.merged.bam \
${each} ${entry}.PTM.SE.sorted.bam \
${entry}.PTT.SE.sorted.bam
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
nucmer -mum -mincluster 65 -minmatch 40 --prefix=${nucmer}${entry} /home/yuzon/references/hybrid_references/PTT_0-1_assembly.v14.fa ${entry}.PTMPTT.fasta
delta-filter -r -q ${nucmer}${entry}.delta >${nucmer}${entry}.filter
show-snps -Clr ${nucmer}${entry}.filter > ${nucmer}${entry}.snps
show-coords ${nucmer}${entry}.delta  > ${nucmer}${entry}.show-coords
grep -v "contig.*PTM" ${nucmer}${entry}.show-coords \
|awk '{print $12 "\t" $13 "\t" $1 "\t" $2 "\t" $4 "\t" $5}'>${nucmer}${entry}.show-coords.txt

nucmer -mum -mincluster 65 -minmatch 40 --prefix=${nucmer}${entry}.alt  /home/yuzon/references/Wyatt_Ellwood_assemblies/15A.fasta ${entry}.PTMPTT.fasta
delta-filter -r -q ${nucmer}${entry}.alt.delta >${nucmer}${entry}.alt.filter
show-snps -Clr ${nucmer}${entry}.alt.filter > ${nucmer}${entry}.alt.snps

cd ${nucmer}

echo "##### Remove Sequencing Error 0-1.PTM.bed and 15.PTT.bed #####"
grep "VBVL" 0-1.show-coords|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai >$entry.0-1_PTM.error
grep -v "VBVL" 15a.show-coords|grep "chrPTT\|contig"|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai >$entry.15a_PTT.error

grep "VBVL" ${entry}.show-coords|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai | bedtools subtract -a stdin -b $entry.0-1_PTM.error > ${entry}.PTM.bed
grep -v "VBVL" ${entry}.show-coords|grep "chrPTT\|contig"|awk '{print $12 "\t" $1 "\t" $2 "\t" $13 ":" $4 "-" $5}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai | bedtools subtract -a stdin -b $entry.15a_PTT.error >${entry}.PTT.bed

echo "######################### regions unique to PTT or PTM ###########################"

bedtools intersect -a ${entry}.PTT.bed -b ${entry}.PTM.bed -nonamecheck \
|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "common"}'> ${entry}.common.bed
bedtools subtract -a ${entry}.PTM.bed -b ${entry}.PTT.bed -nonamecheck \
|cut -f1,2,3|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "PTM"}'> ${entry}.PTMuniq.bed
bedtools subtract -a ${entry}.PTT.bed -b ${entry}.PTM.bed -nonamecheck \
|bedtools sort -i stdin -faidx ${refPTT}.fai|bedtools merge -i stdin|awk   '{print $0 "\t" "PTT"}'> ${entry}.PTTuniq.bed

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
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTM"'|awk -v OFS='\t' '{print $0, "PTM"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTM"}' >$entry.ptm.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "PTT"'|awk -v OFS='\t' '{print $0, "PTT"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "PTT"}' >$entry.ptt.tmp
awk -v OFS='\t' '{print $1,$2,$3,$4}' $entry.tmp| awk  -v OFS='\t' '$4 == "common"'|awk -v OFS='\t' '{print $0, "common"}' \
|bedtools merge -i stdin -d 1 |awk -v OFS='\t' '{print $0, "common"}' >$entry.common.tmp

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
awk -v OFS='\t' '{print $14,$1,$1,$15}' ${entry}.snps |grep "V" |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.PTM.snps
awk -v OFS='\t' '{print $15,$4,$4,$15}' ${entry}.alt.snps |grep "chr\|contig" |bedtools sort -i stdin -faidx ${refPTT}.fai > ${entry}.PTT.snps
bedtools subtract -a ${entry}.PTM.snps -b ${entry}.PTT.snps -nonamecheck |bedtools intersect -a $entry.ptm.tmp -b stdin -wa -nonamecheck \
|bedtools groupby -i stdin -c 1 -o count -nonamecheck |awk '$4 > 5' | awk -v variable="$entry" '{print $1 "\t" $2 "\t" $3 "\t" "PTM" "\t" variable "\t" $4}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai  >${entry}.PTM.support
bedtools subtract -a ${entry}.PTT.snps -b ${entry}.PTM.snps -nonamecheck |bedtools intersect -a $entry.ptt.tmp -b stdin -wa -nonamecheck \
|bedtools groupby -i stdin -c 1 -o count -nonamecheck |awk '$4 > 5' | awk -v variable="$entry" '{print $1 "\t" $2 "\t" $3 "\t" "PTT" "\t" variable "\t" $4}' \
|bedtools sort -i stdin -faidx ${refPTT}.fai  >${entry}.PTT.support


cat ${entry}.PTM.support ${entry}.PTT.support $entry.common.tmp|cut -f1-4|bedtools sort -i stdin -faidx ${refPTT}.fai \
| awk '{print $0 "\t" $3-$2}' | awk -v OFS="\t" '$5>=10000{print $1,$2,$3,$4}'> ${entry}.common_unique.blocks.bed

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
	     else if($1==$5 && $4=="PTM" && $8=="PTT"){print $1, $6, $2, "break"} \
	     else if($1==$5 && $4=="PTT" && $8=="PTM"){print $1, $6, $2, "break"}}' ${entry}.updown.tmp \
| awk -v variable="$entry" '{print $0 "\t" variable} '  >$entry.break
