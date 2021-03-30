#!/bin/bash

#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=PTMxPTMjoint
#  how many cpus are requested
#SBATCH --ntasks=4
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=70G
#  write std out and std error to these files
#SBATCH --error=PTMxPTMjoint.%J.err
#SBATCH --output=PTMxPTMjoint.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

########################################################################
refPTT=/home/yuzon/references/hybrid_references/PTM_FGOB10Ptm-1_assembly.v7.fasta
refPTM=/home/yuzon/references/Wyatt_Ellwood_assemblies/SG1.fasta
snps=/home/yuzon/FGOB10Ptm-1_X_SG1_ptm/SNPcalling/
parsnp_qry=/home/yuzon/FGOB10Ptm-1_X_SG1_ptm/SNPcalling/parsnp_qry/
parsnp_out=/home/yuzon/FGOB10Ptm-1_X_SG1_ptm/SNPcalling/parsnp_out/

#run bash to submit parallel jobs on the cluster
### example: sbatch SL_PTTxPTMhybrid.parallel.sh 100 /home/yuzon/PTTxPTMtest/

########################################################################
cd SNPcalling
cat \
*.PTMPTT.flt.bed \
>temp
sed 's/mitochondrial_contig/PTM	mitochondria/g' temp \
|sed 's/0-1_contig_m86/PTT	mitochondria/g' \
|sed 's/0-1_contig_/PTT	contig/g' \
|sed 's/_/        chr/g' |sed 's/chrPTM/PTM/g' \
|sed 's/chrPTT/PTT/g' |grep "+" \
|awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' \
>all.PTMPTT.flt.bed
grep -v "contig" all.PTMPTT.flt.bed >all.PTMPTT.flt.nuclmito.bed

bedtools intersect -a temp -b matloci.bed -nonamecheck> test
grep 'mito' all.PTMPTT.flt.nuclmito.bed|cat test - >all.matmito.bed

##########Intermediate files#################
cd nucmer
for entry in *common_uniq.bed
do
	awk -v variable="$entry" '{print $0 "\t" variable} ' $entry| sed 's/.common_uniq.bed//g' >$entry.tmp
done

cat *common_uniq.bed.tmp|grep -v "all" >all.common_uniq.bed
bedtools intersect -a all.common_uniq.bed -b matloci.bed -nonamecheck > test
grep 'mito\|m86' all.common_uniq.bed|cat test - >all.common_uniq.matmito.bed

sed 's/PTT_//g' all.common_uniq.bed| sed 's/0-1_contig_m86/mitochondria/g' >test
mv test all.common_uniq.bed
sed 's/PTT_//g' all.common_uniq.matmito.bed| sed 's/0-1_contig_m86/mitochondria/g' >test
mv test all.common_uniq.matmito.bed

############Block files#####################
for entry in *common_unique.blocks.bed
do
  	awk -v variable="$entry" '{print $0 "\t" variable} ' $entry| sed 's/.common_unique.blocks.bed//g' >$entry.tmp
done


cat \
FxS1.common_unique.blocks.bed.tmp \
FxS2.common_unique.blocks.bed.tmp \
FxS3.common_unique.blocks.bed.tmp \
FxS4.common_unique.blocks.bed.tmp \
FxS5.common_unique.blocks.bed.tmp \
FxS6.common_unique.blocks.bed.tmp \
FxS7.common_unique.blocks.bed.tmp \
FxS8.common_unique.blocks.bed.tmp \
FxS9.common_unique.blocks.bed.tmp \
FxS10.common_unique.blocks.bed.tmp \
FxS11.common_unique.blocks.bed.tmp \
FxS12.common_unique.blocks.bed.tmp \
FxS13.common_unique.blocks.bed.tmp \
FxS14.common_unique.blocks.bed.tmp \
FxS15.common_unique.blocks.bed.tmp \
FxS16.common_unique.blocks.bed.tmp \
FxS17.common_unique.blocks.bed.tmp \
FxS18.common_unique.blocks.bed.tmp \
FxS19.common_unique.blocks.bed.tmp \
FxS20.common_unique.blocks.bed.tmp \
FxS21.common_unique.blocks.bed.tmp \
FxS22.common_unique.blocks.bed.tmp \
FxS24.common_unique.blocks.bed.tmp \
FxS25.common_unique.blocks.bed.tmp \
FxS26.common_unique.blocks.bed.tmp \
FxS27.common_unique.blocks.bed.tmp \
FxS28.common_unique.blocks.bed.tmp \
FxS29.common_unique.blocks.bed.tmp \
FxS31.common_unique.blocks.bed.tmp \
FxS32.common_unique.blocks.bed.tmp \
FxS33.common_unique.blocks.bed.tmp \
FxS34.common_unique.blocks.bed.tmp \
FxS35.common_unique.blocks.bed.tmp \
FxS36.common_unique.blocks.bed.tmp \
FxS37.common_unique.blocks.bed.tmp \
FxS39.common_unique.blocks.bed.tmp \
FxS40.common_unique.blocks.bed.tmp \
FxS41.common_unique.blocks.bed.tmp \
FxS42.common_unique.blocks.bed.tmp \
FxS43.common_unique.blocks.bed.tmp \
FxS44.common_unique.blocks.bed.tmp \
FxS45.common_unique.blocks.bed.tmp \
FxS46.common_unique.blocks.bed.tmp \
FxS47.common_unique.blocks.bed.tmp \
FxS48.common_unique.blocks.bed.tmp \
FxS49.common_unique.blocks.bed.tmp \
FxS50.common_unique.blocks.bed.tmp \
FxS51.common_unique.blocks.bed.tmp \
FxS52.common_unique.blocks.bed.tmp \
FxS53.common_unique.blocks.bed.tmp \
FxS54.common_unique.blocks.bed.tmp \
FxS55.common_unique.blocks.bed.tmp \
FxS56.common_unique.blocks.bed.tmp \
FxS57.common_unique.blocks.bed.tmp \
FxS58.common_unique.blocks.bed.tmp \
FxS59.common_unique.blocks.bed.tmp \
FxS60.common_unique.blocks.bed.tmp \
FxS61.common_unique.blocks.bed.tmp \
FxS62.common_unique.blocks.bed.tmp \
FxS66.common_unique.blocks.bed.tmp \
FxS67.common_unique.blocks.bed.tmp \
FxS68.common_unique.blocks.bed.tmp \
FxS69.common_unique.blocks.bed.tmp \
FxS70.common_unique.blocks.bed.tmp \
FxS71.common_unique.blocks.bed.tmp \
FxS72.common_unique.blocks.bed.tmp \
FxS74.common_unique.blocks.bed.tmp \
FxS75.common_unique.blocks.bed.tmp \
FxS76.common_unique.blocks.bed.tmp \
FxS77.common_unique.blocks.bed.tmp \
FxS78.common_unique.blocks.bed.tmp \
FxS80.common_unique.blocks.bed.tmp \
FxS82.common_unique.blocks.bed.tmp \
FxS83.common_unique.blocks.bed.tmp \
FxS84.common_unique.blocks.bed.tmp \
FxS85.common_unique.blocks.bed.tmp \
FxS86.common_unique.blocks.bed.tmp \
FxS87.common_unique.blocks.bed.tmp \
FxS88.common_unique.blocks.bed.tmp \
FxS89.common_unique.blocks.bed.tmp \
FxS90.common_unique.blocks.bed.tmp \
FxS91.common_unique.blocks.bed.tmp \
FxS92.common_unique.blocks.bed.tmp \
FxS93.common_unique.blocks.bed.tmp \
FxS94.common_unique.blocks.bed.tmp \
FxS95.common_unique.blocks.bed.tmp \
FxS96.common_unique.blocks.bed.tmp \
FxS97.common_unique.blocks.bed.tmp \
FxS98.common_unique.blocks.bed.tmp \
FxS99.common_unique.blocks.bed.tmp \
FxS100.common_unique.blocks.bed.tmp \
FxS101.common_unique.blocks.bed.tmp \
FxS102.common_unique.blocks.bed.tmp \
FxS103.common_unique.blocks.bed.tmp \
FxS104.common_unique.blocks.bed.tmp \
FxS105.common_unique.blocks.bed.tmp \
FxS106.common_unique.blocks.bed.tmp \
FxS107.common_unique.blocks.bed.tmp \
FxS108.common_unique.blocks.bed.tmp \
FxS109.common_unique.blocks.bed.tmp \
FxS110.common_unique.blocks.bed.tmp \
FxS111.common_unique.blocks.bed.tmp \
FxS112.common_unique.blocks.bed.tmp \
FxS113.common_unique.blocks.bed.tmp \
FxS114.common_unique.blocks.bed.tmp \
FxS115.common_unique.blocks.bed.tmp \
FxS116.common_unique.blocks.bed.tmp \
FxS117.common_unique.blocks.bed.tmp \
FxS118.common_unique.blocks.bed.tmp \
|grep -v "all" >all.common_unique.blocks.bed
bedtools intersect -a all.common_unique.blocks.bed -b matloci.bed -nonamecheck> test
grep 'mito\|m86' all.common_unique.blocks.bed|cat test - >all.common_unique.blocks.matmito.bed

sed "s/PTM_//g" all.common_unique.blocks.bed >all.common_unique.blocks.bed.tmp

mkdir Rplotfiles
cp all.common* Rplotfiles
cp ../all.* Rplotfiles
################# Segregation Distortion Files #####################
##awk '$6==1 & $7==1 {print $0}' SNPcalling/nucmer/all.recomblocks.bed ##check if error such that PTT and PTM alles in an isolate overlap (begin with even=PTM alleles odd=PTT alleles)
cut -f1,2 all.common_unique.blocks.bed| sort -k1,1 -k2,2n  >start
cut -f1,3 all.common_unique.blocks.bed| sort -k1,1 -k2,2n  >stop
cat start stop | sort -k1,1 -k2,2n |uniq>tmp.ranges
tail -n +2 tmp.ranges|cut -f1,2|paste tmp.ranges -|awk '$1==$3 {print $1 "\t" $2 "\t" $4}'|bedtools sort -i stdin -faidx $refPTT.fai >ranges.bed

for each in *PTM.support
do
entry=${each%.PTM.support}
bedtools intersect -a ranges.bed -b $entry.common_unique.blocks.bed -wao -nonamecheck > $entry.PTM.ranges.bed
bedtools intersect -a ranges.bed -b $entry.common_unique.blocks.bed -wao -nonamecheck > $entry.PTT.ranges.bed
echo -e "chrom\tstart\t stop\t$entry.PTM.support\t$entry.PTT.support" >header
paste  $entry.PTM.ranges.bed  $entry.PTT.ranges.bed |cut -f1-3,7,15 | cat header - |cut -f4,5> $entry.ranges
#sed 's/PTT_/ptt_/g' $entry.PTM.ranges.bed|cat $entry.PTT.ranges.bed -|cut -f7| (echo $entry && cat -) > $entry.merge.ranges

awk '{if ($1 ~ /support/) {print $1;} \
else if($1=="." && $2=="."){ print ".";} \
else if($1=="PTM" && $2=="."){ print "PTM";} \
else if($1=="." && $2=="PTT"){ print "PTT";} \
else if($1=="PTM" && $2=="PTT"){ print ".";}}' $entry.ranges \
|sed 's/.PTM.support//g' > $entry.merge.ranges
done

#Remove all “failed” sequences --less than 50% coverage
cut -f1-3 header|cat - ranges.bed| paste - \
FxS1.ranges \
FxS2.ranges \
FxS3.ranges \
FxS4.ranges \
FxS5.ranges \
FxS6.ranges \
FxS7.ranges \
FxS8.ranges \
FxS9.ranges \
FxS10.ranges \
FxS11.ranges \
FxS12.ranges \
FxS13.ranges \
FxS14.ranges \
FxS15.ranges \
FxS16.ranges \
FxS17.ranges \
FxS18.ranges \
FxS19.ranges \
FxS20.ranges \
FxS21.ranges \
FxS22.ranges \
FxS24.ranges \
FxS25.ranges \
FxS26.ranges \
FxS27.ranges \
FxS28.ranges \
FxS29.ranges \
FxS31.ranges \
FxS32.ranges \
FxS33.ranges \
FxS34.ranges \
FxS35.ranges \
FxS36.ranges \
FxS37.ranges \
FxS39.ranges \
FxS40.ranges \
FxS41.ranges \
FxS42.ranges \
FxS43.ranges \
FxS44.ranges \
FxS45.ranges \
FxS46.ranges \
FxS47.ranges \
FxS48.ranges \
FxS49.ranges \
FxS50.ranges \
FxS51.ranges \
FxS52.ranges \
FxS53.ranges \
FxS54.ranges \
FxS55.ranges \
FxS56.ranges \
FxS57.ranges \
FxS58.ranges \
FxS59.ranges \
FxS60.ranges \
FxS61.ranges \
FxS62.ranges \
FxS66.ranges \
FxS67.ranges \
FxS68.ranges \
FxS69.ranges \
FxS70.ranges \
FxS71.ranges \
FxS72.ranges \
FxS74.ranges \
FxS75.ranges \
FxS76.ranges \
FxS77.ranges \
FxS78.ranges \
FxS80.ranges \
FxS82.ranges \
FxS83.ranges \
FxS84.ranges \
FxS85.ranges \
FxS86.ranges \
FxS87.ranges \
FxS88.ranges \
FxS89.ranges \
FxS90.ranges \
FxS91.ranges \
FxS92.ranges \
FxS93.ranges \
FxS94.ranges \
FxS95.ranges \
FxS96.ranges \
FxS97.ranges \
FxS98.ranges \
FxS99.ranges \
FxS100.ranges \
FxS101.ranges \
FxS102.ranges \
FxS103.ranges \
FxS104.ranges \
FxS105.ranges \
FxS106.ranges \
FxS107.ranges \
FxS108.ranges \
FxS109.ranges \
FxS110.ranges \
FxS111.ranges \
FxS112.ranges \
FxS113.ranges \
FxS114.ranges \
FxS115.ranges \
FxS116.ranges \
FxS117.ranges \
FxS118.ranges \
>all.recomblocks.bed
sed 's/PTM_//g' all.recomblocks.bed>all.recomblocks.bed.tmp

paste \
1.merge.ranges \
2.merge.ranges \
3.merge.ranges \
5.merge.ranges \
6.merge.ranges \
7.merge.ranges \
8.merge.ranges \
9.merge.ranges \
10.merge.ranges \
11.merge.ranges \
12.merge.ranges \
13.merge.ranges \
14.merge.ranges \
15.merge.ranges \
16.merge.ranges \
17.merge.ranges \
18.merge.ranges \
19.merge.ranges \
20.merge.ranges \
21.merge.ranges \
22.merge.ranges \
23.merge.ranges \
25.merge.ranges \
27.merge.ranges \
28.merge.ranges \
29.merge.ranges \
30.merge.ranges \
31.merge.ranges \
32.merge.ranges \
33.merge.ranges \
34.merge.ranges \
36.merge.ranges \
37.merge.ranges \
38.merge.ranges \
39.merge.ranges \
40.merge.ranges \
41.merge.ranges \
42.merge.ranges \
43.merge.ranges \
45.merge.ranges \
46.merge.ranges \
47.merge.ranges \
48.merge.ranges \
49.merge.ranges \
50.merge.ranges \
51.merge.ranges \
52.merge.ranges \
53.merge.ranges \
55.merge.ranges \
56.merge.ranges \
57.merge.ranges \
58.merge.ranges \
59.merge.ranges \
60.merge.ranges \
61.merge.ranges \
62.merge.ranges \
63.merge.ranges \
64.merge.ranges \
65.merge.ranges \
66.merge.ranges \
67.merge.ranges \
68.merge.ranges \
69.merge.ranges \
70.merge.ranges \
72.merge.ranges \
75.merge.ranges \
76.merge.ranges \
77.merge.ranges \
79b.merge.ranges \
81.merge.ranges \
82.merge.ranges \
83.merge.ranges \
84.merge.ranges \
85.merge.ranges \
86.merge.ranges \
87.merge.ranges \
88.merge.ranges \
90.merge.ranges \
91.merge.ranges \
92.merge.ranges \
93.merge.ranges \
95.merge.ranges \
96.merge.ranges \
97.merge.ranges \
98.merge.ranges \
100.merge.ranges \
101.merge.ranges \
102.merge.ranges \
103.merge.ranges \
104.merge.ranges \
105.merge.ranges \
106.merge.ranges \
107.merge.ranges \
108.merge.ranges \
109.merge.ranges \
110.merge.ranges \
111.merge.ranges \
112.merge.ranges \
113.merge.ranges \
114.merge.ranges \
115.merge.ranges \
116.merge.ranges \
117.merge.ranges \
118.merge.ranges \
119.merge.ranges \
120.merge.ranges \
>all.recomblocks.vcf.tmp
rm *merge.ranges
####################################################### LD: r2 and r ##################################################################
### convert file format to vcf###
grep "##" ../32.PTMPTT.flt.vcf >header
awk -v OFS="\t" '{print $1,$2,$1":"$2"-"$3,"T","A",".",".",".","GT:AD:DP:GQ:PL"}' all.recomblocks.bed >ptt.columns
sed 's/PTT_/ptt_/g' ptt.columns| sed 's/0-1_contig/ptm/g'|tail -n +2 >ptm.columns
cat ptt.columns ptm.columns \
|sed  "1s/.*/#CHROM POS ID	 REF ALT QUAL FILTER INFO FORMAT/" \
|awk '{$1=$1}1' OFS="\t" >columns
sed 's/PTT/0\/0/g' all.recomblocks.vcf.tmp |sed 's/PTM/0\/1/g' | sed 's/\./\.\/\./g'|paste columns - | cat header - >all.recomblocks.ptm.vcf
sed 's/PTT/0\/1/g' all.recomblocks.vcf.tmp |sed 's/PTM/0\/0/g' | sed 's/\./\.\/\./g'|paste columns - | cat header - >all.recomblocks.ptt.vcf

vcftools --vcf all.recomblocks.ptm.vcf --max-missing-count 100 --out all.recomblocks.ptm --recode
vcftools --vcf all.recomblocks.ptt.vcf --max-missing-count 100 --out all.recomblocks.ptt --recode
#####check if error such that PTT and PTM alles in an isolate overlap
grep 'chrPTT_.*chrPTT_.*PTM.*chrPTT_.*PTT' *.ranges.bed
grep 'contig.*contig.*PTM.*contig.*PTT' *.ranges.bed

#Create a r2 and r matrix 
##Plink
/data/biosoftware/plink/v1.90beta4/plink --r square --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r.ptm.matrix --allow-extra-chr
/data/biosoftware/plink/v1.90beta4/plink --r2 square --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r2.ptm.matrix --allow-extra-chr

/data/biosoftware/plink/v1.90beta4/plink --r square --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r.ptt.matrix --allow-extra-chr
/data/biosoftware/plink/v1.90beta4/plink --r2 square --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r2.ptt.matrix --allow-extra-chr

/data/biosoftware/plink/v1.90beta4/plink --r2 inter-chr --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r2.ptm --allow-extra-chr
/data/biosoftware/plink/v1.90beta4/plink --r inter-chr --vcf all.recomblocks.ptm.recode.vcf --recode structure --out all.r.ptm --allow-extra-chr

/data/biosoftware/plink/v1.90beta4/plink --r2 inter-chr --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r2.ptt --allow-extra-chr
/data/biosoftware/plink/v1.90beta4/plink --r inter-chr --vcf all.recomblocks.ptt.recode.vcf --recode structure --out all.r.ptt --allow-extra-chr

awk '$1!=$4  {print $0}' all.r.ptm.ld >all.r.ptm.inter.ld
awk '$1!=$4  {print $0}' all.r2.ptm.ld >all.r2.ptm.inter.ld

awk '$1!=$4  {print $0}' all.r.ptt.ld >all.r.ptt.inter.ld
awk '$1!=$4  {print $0}' all.r2.ptt.ld >all.r2.ptt.inter.ld

module load R
Rscript --vanilla LinkageDisequilibrium.R all.r.ptm.matrix.ld all.r2.ptm.matrix.ld all.r2.ptm.ld all.r.ptm.inter.ld
Rscript --vanilla LinkageDisequilibrium.R all.r.ptt.matrix.ld all.r2.ptt.matrix.ld all.r2.ptt.ld all.r.ptt.inter.ld


####################################################### Recombination Breakpoints ##################################################################
cat \
FxS1.break \
FxS2.break \
FxS3.break \
FxS4.break \
FxS5.break \
FxS6.break \
FxS7.break \
FxS8.break \
FxS9.break \
FxS10.break \
FxS11.break \
FxS12.break \
FxS13.break \
FxS14.break \
FxS15.break \
FxS16.break \
FxS17.break \
FxS18.break \
FxS19.break \
FxS20.break \
FxS21.break \
FxS22.break \
FxS24.break \
FxS25.break \
FxS26.break \
FxS27.break \
FxS28.break \
FxS29.break \
FxS31.break \
FxS32.break \
FxS33.break \
FxS34.break \
FxS35.break \
FxS36.break \
FxS37.break \
FxS39.break \
FxS40.break \
FxS41.break \
FxS42.break \
FxS43.break \
FxS44.break \
FxS45.break \
FxS46.break \
FxS47.break \
FxS48.break \
FxS49.break \
FxS50.break \
FxS51.break \
FxS52.break \
FxS53.break \
FxS54.break \
FxS55.break \
FxS56.break \
FxS57.break \
FxS58.break \
FxS59.break \
FxS60.break \
FxS61.break \
FxS62.break \
FxS66.break \
FxS67.break \
FxS68.break \
FxS69.break \
FxS70.break \
FxS71.break \
FxS72.break \
FxS74.break \
FxS75.break \
FxS76.break \
FxS77.break \
FxS78.break \
FxS80.break \
FxS82.break \
FxS83.break \
FxS84.break \
FxS85.break \
FxS86.break \
FxS87.break \
FxS88.break \
FxS89.break \
FxS90.break \
FxS91.break \
FxS92.break \
FxS93.break \
FxS94.break \
FxS95.break \
FxS96.break \
FxS97.break \
FxS98.break \
FxS99.break \
FxS100.break \
FxS101.break \
FxS102.break \
FxS103.break \
FxS104.break \
FxS105.break \
FxS106.break \
FxS107.break \
FxS108.break \
FxS109.break \
FxS110.break \
FxS111.break \
FxS112.break \
FxS113.break \
FxS114.break \
FxS115.break \
FxS116.break \
FxS117.break \
FxS118.break \
|bedtools sort -i stdin -faidx $refPTT.fai \
>all.breaks

bedtools makewindows -g $refPTT.fai -w 10000 -s 10000 > genome.windows.bed
bedtools intersect -a genome.windows.bed -b all.breaks -c -nonamecheck \
|bedtools intersect -nonamecheck -wa -wb -a stdin -b all.breaks \
|awk '{print $0 "\t" $7-$6}'  > all.breaks.density

####### Get Parental and Hybrid Frequencies
vcftools --vcf all.recomblocks.ptm.recode.vcf --out all.recomblocks.ptm.recode --plink

while read p
do
  pairs=`echo "$p"| grep -v "SNP"|awk '{print $3 "\t" $6}'`
  /data/biosoftware/plink/v1.07/plink-1.07-x86_64/plink --ld $pairs --file all.recomblocks.ptm.recode --noweb \
  |grep "AA\|TT\|AT\|TA\|A0\|T0\|0A\|0T" |awk '{print $2}' |grep -v "In\|LD\|ld\|info\|phase"| awk 'BEGIN { ORS = " " } { print }'
  echo "$pairs"
done <all.r.ptm.ld >all.pairs.freq



