#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=Populations
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=120:00:00
#  maximum requested memory
#SBATCH --mem=75G
#  write std out and std error to these files
#SBATCH --error=Populations.%J.err
#SBATCH --output=Populations.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=global

###########################################################
echo "Software"
###########################################################
##module load java/x64/8u121

merged_genome=/home/yuzon/references/hybrid_references/merged_genome_PTTxPTM.fasta

harvest=/data/biosoftware/harvest/harvest
Sorted=/home/yuzon/Populations/aligned/sorted/


##########################################################
echo "Remove - in squence names"
##########################################################
for each in /home/yuzon/Populations/SNPcalling/parsnp_qry/mitochondria/*fasta
do
	sed -i 's/-//g' $each
done

##########################################################
##### SNP call for Illumina Data
cd $Sorted

samtools mpileup -q 30 -Q15 -uf ${merged_genome} -g -t DP,AD,ADF,ADR,SP --skip-indels -L 999999999 -d 999999999 \
E32_Pt54m3.PTMPTT.read.bam \
E74S_Pt63_AB.PTMPTT.read.bam \
G106_Pt124m3.PTMPTT.read.bam \
H165_Pt152m3.PTMPTT.read.bam \
M207_Pt97m2.PTMPTT.read.bam \
T305_Pt88m2.PTMPTT.read.bam \
E47S_Pt39.PTMPTT.read.bam \
T81_Pt74.PTMPTT.read.bam \
M209_Pt99.PTMPTT.read.bam \
M211_Pt101.PTMPTT.read.bam \
M218_Pt103.PTMPTT.read.bam \
M219_Pt104.PTMPTT.read.bam \
M95_Pt119_AB.PTMPTT.read.bam \
M221_Pt105.PTMPTT.read.bam \
E64S_Pt46.PTMPTT.read.bam \
M208_Pt98.PTMPTT.read.bam \
M210_Pt100.PTMPTT.read.bam \
M97_Pt118.PTMPTT.read.bam \
G11S_Pt131.PTMPTT.read.bam \
G101_Pt145.PTMPTT.read.bam \
H479_Pt171.PTMPTT.read.bam \
H472_Pt169_AB.PTMPTT.read.bam \
H470_Pt177_AB.PTMPTT.read.bam \
Pgraminea_Illumina.read.bam \
|bcftools call -m -Ov --ploidy 1 --skip-variants indels -o /home/yuzon/Populations/SNPcalling/Illumina.vcf

cd /home/yuzon/Populations/SNPcalling/

bgzip -c -f Illumina.vcf >Illumina.vcf.gz
tabix -p vcf -f Illumina.vcf.gz

#### Filter VCF file
bcftools filter -e 'FORMAT/AD[0] < 5 & FORMAT/GT = "0" | FORMAT/AD[1] < 5 & FORMAT/GT = "1" | FORMAT/AD[2] < 5 & FORMAT/GT = "2"'  --set-GTs .  Illumina.vcf.gz \
| bcftools filter -e ' ((FORMAT/AD[0])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "0" |  ((FORMAT/AD[1])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "1"|  ((FORMAT/AD[2])/(FORMAT/DP)) < 0.7 & FORMAT/GT = "2"'  --set-GTs . \
| bcftools filter -e ' %QUAL<20 || MQ < 30 '  >Illumina.flt.vcf

bgzip -f -c Illumina.flt.vcf >Illumina.flt.vcf.gz
tabix -f -p vcf Illumina.flt.vcf.gz

#### SNP call for Genome Assembly Data
${harvest}/parsnp -c -v -r /home/yuzon/references/parsnp_ref/test.fa -d /home/yuzon/Populations/SNPcalling/parsnp_qry/parsnp_qry_2/ -o  /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2
${harvest}/harvesttools -i /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.ggr -V /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.vcf
grep "#\|PASS" /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.vcf >/home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.PASS.vcf
vcftools --vcf /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.vcf --extract-FORMAT-info GT --out /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.nucl

echo "##fileformat=VCFv4.2"|cat - /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.vcf \
|grep '#\|PASS' \
|sed 's/NA/\./g' \
|sed 's/chr/chrPTT_/g' \
|sed 's/01_contig/0-1_contig/g' \
 >/home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.format.vcf

bgzip -f -c /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.format.vcf > /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.format.vcf.gz
tabix -f -p vcf /home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.format.vcf.gz

/data/biosoftware/bcftools/bcftools-1.4/bcftools merge \
/home/yuzon/Populations/SNPcalling/parsnp_out/parsnp_out_2/parsnp.format.vcf.gz \
/home/yuzon/Populations/SNPcalling/Illumina.flt.vcf.gz \
-Oz -o /home/yuzon/Populations/SNPcalling/phylo_map_to_PTT/Merged.vcf.gz

### keep loci with at least one snp allele
cd phylo_map_to_PTT

/data/biosoftware/bcftools/bcftools-1.4/bcftools view --types snps \
/home/yuzon/Populations/SNPcalling/phylo_map_to_PTT/Merged.vcf.gz \
-Oz -o /home/yuzon/Populations/SNPcalling/phylo_map_to_PTT/Merged.snps.vcf.gz
tabix -f /home/yuzon/Populations/SNPcalling/phylo_map_to_PTT/Merged.snps.vcf.gz

vcftools --max-missing 1.00 --recode  \
--gzvcf Merged.snps.vcf.gz \
--out Merged.mis

bcftools view -s ^T305_Pt88m2.PTMPTT Merged.mis.recode.vcf \
|/data/biosoftware/bcftools/bcftools-1.4/bcftools view --types snps -c1 >tmp.vcf
mv tmp.vcf Merged.mis.recode.vcf

###remove loci that have non-SNP alleles
awk '$5 !~ /AA|AC|AG|AT|CA|CC|CG|CT|GA|GC|GG|GT|TA|TC|TG|TT/ {print $0}' Merged.mis.recode.vcf|grep -v 'LowQual.*GT:' > Merged.flt.recode.vcf 

bgzip -c -f Merged.flt.recode.vcf >Merged.flt.recode.vcf.gz
tabix -p -f vcf Merged.flt.recode.vcf.gz

### missingness (GT) per individual
vcftools --gzvcf Merged.snps.vcf.gz \
--missing-indv \
--out Merged.mis.indv  #most missing GT: T305_Pt88m2 (71.5%)

#### Format VCF to FASTA
module load perl
zcat Merged.flt.recode.vcf.gz |/data/biosoftware/vcftools/vcftools-0.1.14/src/perl/vcf-to-tab \
|sed 's/\///g' >Merged.flt.GT.FORMAT
cut -f4-  Merged.flt.GT.FORMAT \
|datamash transpose |awk '{print ">" $0}' \
| sed 's/test.fa.ref/0-1\n/g'| sed 's/\.PTMPTT/\n/g'|sed 's/\.fasta/\n/g' \
|sed 's/\.bam/\n/g' | sed 's/\./-/g' |sed 's/\-fna/\n/g' \
|sed 's/Pyrenophora_tritici-repentis_Pt/>Ptr/g' \
|sed 's/.*_Ptr/>Ptr/g' |sed 's/.*_PTR/>Ptr/g' |sed 's/_genomic//g' \
|sed 's/Pgraminea_Illumina-sorted/>Pg_QWC/g' \
|sed 's/PTM_FGOB10Ptm-1_assembly-v7/>FGOB10Ptm-1/g' \
|tr -d " \t" >Merged.flt.fasta

echo "############ Divergence Stats (Dxy and Fst): only PTT and PTM isolates #########################"
bcftools view -S ptm_samplelist.txt Merged.snps.vcf.gz >Merged.flt.recode.ptm.vcf.gz
bcftools view -S ptt_samplelist.txt Merged.snps.vcf.gz >Merged.flt.recode.ptt.vcf.gz 

vcftools --max-missing 0.7 --recode  --gzvcf  Merged.flt.recode.ptm.vcf.gz --out Merged.ptm.flt2
vcftools --max-missing 0.7 --recode  --gzvcf  Merged.flt.recode.ptt.vcf.gz --out Merged.ptt.flt2

bgzip -f -c Merged.ptm.flt2.recode.vcf >Merged.ptm.flt2.recode.vcf.gz
bgzip -f -c Merged.ptt.flt2.recode.vcf >Merged.ptt.flt2.recode.vcf.gz
tabix -f Merged.ptm.flt2.recode.vcf.gz
tabix -f Merged.ptt.flt2.recode.vcf.gz

bcftools merge Merged.ptm.flt2.recode.vcf.gz Merged.ptt.flt2.recode.vcf.gz -Ov > Merged.flt.dxy.tmp.vcf

/data/biosoftware/bcftools/bcftools-1.4/bcftools view --types snps -c1 \
Merged.flt.dxy.tmp.vcf > Merged.flt.dxy.vcf


echo "############ Get Derived/Polarized SNPs #########################"
#bcftools view -S ptr_samplelist.txt Merged.snps.vcf.gz |bcftools filter -e 'AC<1' \
#|grep -v '	1:.*	0:\|	0:.*	1:\|	2:.*	0:\|	0:.*	2:\|	2:.*	1:\|	1:.*	2:' \
#> Merged.flt.recode.ptr.vcf
#bgzip -f Merged.flt.recode.ptr.vcf>Merged.flt.recode.ptr.vcf.gz
#tabix -f Merged.flt.recode.ptr.vcf.gz
#cut -f1 pops.abbababa.txt>ptrpttpg.txt
#bcftools view --regions-file Merged.flt.recode.ptr.vcf.gz Merged.snps.vcf.gz \
#|bcftools view -m2 -M2 -v snps -S ptrpttpg.txt |sed 's/PTT_//g'>Derived.bi.abbababa.vcf


