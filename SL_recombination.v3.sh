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
#SBATCH --time=80:00:00
#  maximum requested memory
#SBATCH --mem=50M #70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTMrecomb.%J.err
#SBATCH --output=PTTxPTMrecomb.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

##########################################################################################################
################### Get Parental and Hybrid Frequencies ##################################################
##########################################################################################################

#### Output: parental, recomb, recomb, parental: AA TA AT TT
date 

file=$1  #all.r.ptm.inter.ld

echo $file

cat $file |while read p
do

region1=`echo "$p"| grep -v "SNP"|awk '{print $3}'`
region2=`echo "$p"| grep -v "SNP"|awk '{print $6}'`
grep "$region1\|$region2\|CHROM"  all.recomblocks.ptm.recode.vcf  |cut -f3,10-|datamash transpose>${region1}_${region2}.tmp
tt=`awk '$2 ~ "0/0"&& $3 ~  "0/0" {print $0}' ${region1}_${region2}.tmp|wc -l`
tm=`awk '$2 ~ "0/0"&& $3 ~  "0/1" {print $0}' ${region1}_${region2}.tmp|wc -l`
mt=`awk '$2 ~ "0/1"&& $3 ~  "0/0" {print $0}' ${region1}_${region2}.tmp|wc -l`
mm=`awk '$2 ~ "0/1"&& $3 ~  "0/1" {print $0}' ${region1}_${region2}.tmp|wc -l`
echo $region1	$region2	$tt	$tm	$mt	$mm |awk -v OFS='\t' '{print $0,($3*$6)/($4*$5),$3/$6}'
rm ${region1}_${region2}.tmp

done > $file.txt

date
