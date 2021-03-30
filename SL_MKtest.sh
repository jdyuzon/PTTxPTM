#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=mktest
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=100:00:00 #100:00:00
#  maximum requested memory
#SBATCH --mem=15G  #70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTMhybrid.%J.err
#SBATCH --output=PTTxPTMhybrid.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

############## CDS to Amino Acid ###################
#/home/yuzon/BioSoft/gffread/gffread -J -y /home/yuzon/references/hybrid_references/PTT_0-1_annotation.gff.AA.fasta \
#-g /home/yuzon/references/hybrid_references/PTT_0-1_assembly.v14.fa \
#/home/yuzon/references/hybrid_references/PTT_0-1_annotation.gff.sorted.txt
############## tblastn CDS to pseudogenomes ####################
module load python/3.7.1
mkdir MKtest
cd MKtest

####get all samples in a file
cp ~/Populations/SNPcalling/*.PTMPTT.fasta .
cp /home/yuzon/references/hybrid_references/PTT_0-1_annotation.gff.AA.fasta .
cp ~/references/Pyrenophora_sisters/Pyrenophora_graminea_CBS33629_1/Pgraminea_GCA_012365135.1_ASM1236513v1_genomic.fna .

ls *PTMPTT.fasta > subjects.txt
###filter by sequence length
ruby ~/BioSoft/filter_queries.rb /home/yuzon/references/hybrid_references/PTT_0-1_annotation.gff.AA.fasta exons_red.fasta 50
###adds terminal stop codon to AA sequence
sed '/^[[:space:]]*$/d' exons_red.fasta|sed ':a;N;$!ba;s/\n>/\*\n>/g' > tmp
mv tmp exons_red.fasta
####index fasta files
 for i in *PTMPTT.fasta
  do
  	makeblastdb -in ${i} -dbtype nucl
  done

makeblastdb -in PTT_0-1_annotation.gff.AA.fasta -dbtype nucl
#### get orthologs
python3 ~/BioSoft/find_orthologs.JDY.py -t -s 1 --refine exons_red.fasta subjects.txt


########Filter Alignment #######################
mkdir MKtest_filter

for alignment in *nucl.fasta
do

prefix=${alignment%.fasta}
echo $alignment
echo $prefix

#### check if Start/Stop codon
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $alignment \
| awk -F '\t'  '($2 ~ /TAA$|TAG$|TGA$/)' |awk -F '\t'  '($2 ~ /^ATG/)' | tr "\t" "\n" >MKtest_filter/${prefix}.stopstart.fasta

#### remove sequences with 0.1 gaps: https://www.biostars.org/p/434389/ (Mensur Dlakic)
python3 ~/BioSoft/fasta_drop.py MKtest_filter/${prefix}.stopstart.fasta MKtest_filter/${prefix}.gap.fasta 0.1
#### clean the fasta files to remove extra newlines
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' MKtest_filter/${prefix}.gap.fasta \
>tmp && mv tmp MKtest_filter/${prefix}.gap.fasta

done


#### discard alignments with less than two sequences
cd MKtest_filter
pwd

find -name '*.fasta' | xargs  wc -l | awk '{if($1 < 4 && index($2, "fasta")>0 ) print $2}' | xargs rm
##### discard alignments that don't have PTT and PTM, PG individuals
find . -type f -print0 | xargs --null grep -Z -L '>M' | xargs --null rm 
find . -type f -print0 | xargs --null grep -Z -L '>E\|>G\|>H\|>T' | xargs --null rm

########################## Calculate MKtest using BioPython #################
cp ~/BioSoft/codonalignmentJDY.py .
cp ~/BioSoft/MKtest.py .

module load python/3.7.1
python3 MKtest.py

sort MK_results.pttptm.txt| sed 's/_nucl.gap.fasta//g' >tmp1
grep 'mRNA' ~/references/hybrid_references/PTT_0-1_annotation.gff.sorted.txt \
|sed 's/ID=//g' | sed 's/_mrna.*/_mrna/g' | sort -k9 \
|awk '{print $1 "\t" $4 "\t" $5 "\t" $9}'>tmp2
join -t $'\t' -1 1 -2 4 tmp1 tmp2 -o 2.1,2.2,2.3,2.4,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 \
|sed 's/PTT_//g' >MK_results.pttptm.txt

sort MK_results.ptmptt.txt| sed 's/_nucl.gap.fasta//g' >tmp1
grep 'mRNA' ~/references/hybrid_references/PTT_0-1_annotation.gff.sorted.txt \
|sed 's/ID=//g' | sed 's/_mrna.*/_mrna/g' | sort -k9 \
|awk '{print $1 "\t" $4 "\t" $5 "\t" $9}'>tmp2
join -t $'\t' -1 1 -2 4 tmp1 tmp2 -o 2.1,2.2,2.3,2.4,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 \
|sed 's/PTT_//g' >MK_results.ptmptt.txt

bedtools intersect -wa -wb -a MK_results.pttptm.txt -b ~/PTTxPTMtest/PTTxPTMcompare/circlize/all.segregdistort.PTTxPTM.txt \
|awk -v  OFS='\t' '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$26,$27}' >MK_results.pttptm.segregdistort.txt

bedtools intersect -wa -wb -a MK_results.ptmptt.txt -b ~/PTTxPTMtest/PTTxPTMcompare/circlize/all.segregdistort.PTTxPTM.txt \
|awk -v  OFS='\t' '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$26,$27}' >MK_results.ptmptt.segregdistort.txt

