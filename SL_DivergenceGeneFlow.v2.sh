#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=DivergenceGeneFlow
#  how many cpus are requested
#SBATCH --ntasks=5
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=15G #70G
#  write std out and std error to these files
#SBATCH --error=DivergenceGeneFlow.%J.err
#SBATCH --output=DivergenceGeneFlow.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=global

module load python


#######################################################
############## Dxy and Fst: Simon Martin ##############
#######################################################
##sed 's/\./N/g' Merged.flt.GT.FORMAT |grep -v "CH" >tmp
##grep 'CH' Merged.flt.GT.FORMAT|cat - tmp >Merged.flt.GT.FORMAT.geno

######format geno file: only PTT and PTM######################
cd /home/yuzon/Populations/SNPcalling/phylo_map_to_PTT

sed 's/PTT_//g' Merged.flt.dxy.vcf >Merged.ptmptt.mis.recode.vcf

python ~/BioSoft/genomics_general-master/VCF_processing/parseVCF.py \
-i Merged.ptmptt.mis.recode.vcf --ploidy 1 -o Merged.pttptm.geno

python popgenWindows.py -w 50000 -s 5000 -m 10 \
-g Merged.pttptm.geno -o divergence.csv \
-f haplo -T 5 --popsFile ptmptt_pops.txt \
--writeFailedWindows \
-p ptt M207_Pt97m2.PTMPTT,M208_Pt98.PTMPTT,M209_Pt99.PTMPTT,M210_Pt100.PTMPTT,M211_Pt101.PTMPTT,M218_Pt103.PTMPTT,M219_Pt104.PTMPTT,M221_Pt105.PTMPTT,M95_Pt119_AB.PTMPTT,M97_Pt118.PTMPTT \
-p ptm E32_Pt54m3.PTMPTT,E47S_Pt39.PTMPTT,E64S_Pt46.PTMPTT,E74S_Pt63_AB.PTMPTT,G101_Pt145.PTMPTT,G106_Pt124m3.PTMPTT,G11S_Pt131.PTMPTT,H165_Pt152m3.PTMPTT,H470_Pt177_AB.PTMPTT,H472_Pt169_AB.PTMPTT,H479_Pt171.PTMPTT,T81_Pt74.PTMPTT

sed 's/PTT_//g' divergence.csv |grep -v 'contig' >divergence.csv.tmp


#######################################################
######## ABBA-BABA and f-stats: Simon Martin ##########
#######################################################

######format geno file######################
python ~/BioSoft/genomics_general-master/VCF_processing/parseVCF.py -i Derived.bi.abbababa.vcf --ploidy 1 -o Merged.bi.pgptmpttptr.csv

########## Sliding window ABBA-BABA ###############
module load python/3.7.1
python3 ~/BioSoft/genomics_general-master/ABBABABAwindows.py \
-g Merged.bi.pgptmpttptr.csv --windType coordinate \
-o abbababa.csv -P1 pg -P2 ptm -P3 ptt -O ptr \
--popsFile pops.abbababa.txt -w 500000 -m 100 -s 100000 --T 1 --writeFailedWindows -f haplo \
--haploid M207_Pt97m2.PTMPTT,M208_Pt98.PTMPTT,M209_Pt99.PTMPTT,M210_Pt100.PTMPTT,M211_Pt101.PTMPTT,M218_Pt103.PTMPTT,M219_Pt104.PTMPTT,M221_Pt105.PTMPTT,M95_Pt119_AB.PTMPTT,M97_Pt118.PTMPTT,E32_Pt54m3.PTMPTT,E47S_Pt39.PTMPTT,E64S_Pt46.PTMPTT,E74S_Pt63_AB.PTMPTT,G101_Pt145.PTMPTT,G106_Pt124m3.PTMPTT,G11S_Pt131.PTMPTT,H165_Pt152m3.PTMPTT,H470_Pt177_AB.PTMPTT,H472_Pt169_AB.PTMPTT,H479_Pt171.PTMPTT,T81_Pt74.PTMPTT,Pgraminea_Illumina.sorted.bam,GCA_003171515.2_CUR_PTRM4_2.1_genomic.fna,GCA_003171545.1_PtrARCrossB10_genomic.fna,GCA_003231325.1_Ptr134_genomic.fna,GCA_003231345.1_Ptr5213_genomic.fna,GCA_003231355.1_Ptr11137_genomic.fna,GCA_003231365.1_Ptr239_genomic.fna,GCA_003231415.2_CUR_PtrDW5_2.1_genomic.fna,GCA_003231425.1_Ptr86-124_genomic.fna,GCA_008692205.1_Ptr_V0001_genomic.fna,Pyrenophora_tritici-repentis_Pt-1C-BFP.fasta


sed 's/PTT_//g' abbababa.csv |grep -v 'contig' >abbababa.csv.tmp

######### Genome wide ABBA-BABA ###################
python /home/yuzon/BioSoft/genomics_general-master/freq.py -g Derived.bi.abbababa.vcf \
-p pg -p ptm -p ptt -p ptr \
--popsFile pops.abbababa.txt --target derived \
-o derFreq.tsv.gz
