#!/bin/bash
#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=PTTxPTTjoint
#  how many cpus are requested
#SBATCH --ntasks=4
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTTjoint.%J.err
#SBATCH --output=PTTxPTTjoint.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

########################################################################
refPTT=~/references/hybrid_references/PTT_0-1_assembly.v14.fa
refPTM=~/references/hybrid_references/PTM_FGOB10Ptm-1_assembly.v7.fasta
snps=~/15Ax0-1_ptt/SNPcalling/
parsnp_qry=~/15Ax0-1_ptt/SNPcalling/parsnp_qry/
parsnp_out=~/P15Ax0-1_ptt/SNPcalling/parsnp_out/

#run bash to submit parallel jobs on the cluster
### example: sbatch SL_PTTxPTMhybrid.parallel.sh 100 /home/yuzon/PTTxPTMtest/
########################################################################
cd SNPcalling
cat \
*.PTMPTT.flt.bed \
>temp
sed 's/VBVL01000020.1/PTM	mitochondria/g' temp \
|sed 's/0-1_contig_m86/PTT	mitochondria/g' \
|sed 's/0-1_contig_/PTT	contig/g' \
| sed 's/VBVL01000001.1/chrPTM_1/g' \
| sed 's/VBVL01000002.1/chrPTM_2/g' \
| sed 's/VBVL01000003.1/chrPTM_3/g' \
| sed 's/VBVL01000004.1/chrPTM_5/g' \
| sed 's/VBVL01000005.1/chrPTM_10/g' \
| sed 's/VBVL01000006.1/chrPTM_6/g' \
| sed 's/VBVL01000007.1/chrPTM_7/g' \
| sed 's/VBVL01000008.1/chrPTM_4/g' \
| sed 's/VBVL01000009.1/chrPTM_8/g' \
| sed 's/VBVL01000010.1/chrPTM_9/g' \
| sed 's/VBVL01000011.1/chrPTM_11/g' \
| sed 's/VBVL01000012.1/chrPTM_12/g' \
|sed 's/_/        chr/g' |sed 's/chrPTM/PTM/g' \
|sed 's/chrPTT/PTT/g' |grep "+" \
|awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' \
>all.PTMPTT.flt.bed
grep -v "contig\|V" all.PTMPTT.flt.bed >all.PTMPTT.flt.nuclmito.bed

bedtools intersect -a temp -b matloci.bed -nonamecheck> test
grep 'mito' all.PTMPTT.flt.nuclmito.bed|cat test - >all.matmito.bed

echo "##########Intermediate files#################"
cd nucmer
for entry in *common_uniq.bed
do
	awk -v variable="$entry" '{print $0 "\t" variable} ' $entry| sed 's/.common_uniq.bed//g' >$entry.tmp
done

cat \
*common_uniq.bed.tmp \
|grep -v "all" >all.common_uniq.bed
bedtools intersect -a all.common_uniq.bed -b matloci.bed -nonamecheck > test
grep 'mito\|m86' all.common_uniq.bed|cat test - >all.common_uniq.matmito.bed

sed 's/PTT_//g' all.common_uniq.bed| sed 's/0-1_contig_m86/mitochondria/g' >test
mv test all.common_uniq.bed
sed 's/PTT_//g' all.common_uniq.matmito.bed| sed 's/0-1_contig_m86/mitochondria/g' >test
mv test all.common_uniq.matmito.bed

echo "############Block files#####################"
for entry in *common_unique.blocks.bed
do
  	awk -v variable="$entry" '{print $0 "\t" variable} ' $entry| sed 's/.common_unique.blocks.bed//g' >$entry.tmp
done

cat \
0-1.common_unique.blocks.bed.tmp \
2.common_unique.blocks.bed.tmp \
3.common_unique.blocks.bed.tmp \
4.common_unique.blocks.bed.tmp \
5.common_unique.blocks.bed.tmp \
6.common_unique.blocks.bed.tmp \
7.common_unique.blocks.bed.tmp \
9.common_unique.blocks.bed.tmp \
10.common_unique.blocks.bed.tmp \
11.common_unique.blocks.bed.tmp \
12.common_unique.blocks.bed.tmp \
13.common_unique.blocks.bed.tmp \
14.common_unique.blocks.bed.tmp \
15a.common_unique.blocks.bed.tmp \
15.common_unique.blocks.bed.tmp \
16.common_unique.blocks.bed.tmp \
17.common_unique.blocks.bed.tmp \
18.common_unique.blocks.bed.tmp \
20.common_unique.blocks.bed.tmp \
21.common_unique.blocks.bed.tmp \
22.common_unique.blocks.bed.tmp \
23.common_unique.blocks.bed.tmp \
24.common_unique.blocks.bed.tmp \
27.common_unique.blocks.bed.tmp \
28.common_unique.blocks.bed.tmp \
29.common_unique.blocks.bed.tmp \
32.common_unique.blocks.bed.tmp \
34.common_unique.blocks.bed.tmp \
35.common_unique.blocks.bed.tmp \
36.common_unique.blocks.bed.tmp \
37.common_unique.blocks.bed.tmp \
38.common_unique.blocks.bed.tmp \
41.common_unique.blocks.bed.tmp \
42.common_unique.blocks.bed.tmp \
43.common_unique.blocks.bed.tmp \
45.common_unique.blocks.bed.tmp \
46.common_unique.blocks.bed.tmp \
49.common_unique.blocks.bed.tmp \
50.common_unique.blocks.bed.tmp \
51.common_unique.blocks.bed.tmp \
52.common_unique.blocks.bed.tmp \
55.common_unique.blocks.bed.tmp \
65.common_unique.blocks.bed.tmp \
67.common_unique.blocks.bed.tmp \
68.common_unique.blocks.bed.tmp \
70.common_unique.blocks.bed.tmp \
71.common_unique.blocks.bed.tmp \
72.common_unique.blocks.bed.tmp \
73.common_unique.blocks.bed.tmp \
74.common_unique.blocks.bed.tmp \
75.common_unique.blocks.bed.tmp \
76.common_unique.blocks.bed.tmp \
77.common_unique.blocks.bed.tmp \
81.common_unique.blocks.bed.tmp \
85.common_unique.blocks.bed.tmp \
87.common_unique.blocks.bed.tmp \
88.common_unique.blocks.bed.tmp \
89.common_unique.blocks.bed.tmp \
90.common_unique.blocks.bed.tmp \
91.common_unique.blocks.bed.tmp \
92.common_unique.blocks.bed.tmp \
94.common_unique.blocks.bed.tmp \
95.common_unique.blocks.bed.tmp \
96.common_unique.blocks.bed.tmp \
99.common_unique.blocks.bed.tmp \
101.common_unique.blocks.bed.tmp \
102.common_unique.blocks.bed.tmp \
105.common_unique.blocks.bed.tmp \
108.common_unique.blocks.bed.tmp \
201.common_unique.blocks.bed.tmp \
203.common_unique.blocks.bed.tmp \
204.common_unique.blocks.bed.tmp \
205.common_unique.blocks.bed.tmp \
206.common_unique.blocks.bed.tmp \
207.common_unique.blocks.bed.tmp \
208.common_unique.blocks.bed.tmp \
210.common_unique.blocks.bed.tmp \
212.common_unique.blocks.bed.tmp \
214.common_unique.blocks.bed.tmp \
215.common_unique.blocks.bed.tmp \
216.common_unique.blocks.bed.tmp \
217.common_unique.blocks.bed.tmp \
219.common_unique.blocks.bed.tmp \
222.common_unique.blocks.bed.tmp \
226.common_unique.blocks.bed.tmp \
227.common_unique.blocks.bed.tmp \
229.common_unique.blocks.bed.tmp \
230.common_unique.blocks.bed.tmp \
233.common_unique.blocks.bed.tmp \
236.common_unique.blocks.bed.tmp \
238.common_unique.blocks.bed.tmp \
239.common_unique.blocks.bed.tmp \
240.common_unique.blocks.bed.tmp \
241.common_unique.blocks.bed.tmp \
242.common_unique.blocks.bed.tmp \
243.common_unique.blocks.bed.tmp \
244.common_unique.blocks.bed.tmp \
245.common_unique.blocks.bed.tmp \
246.common_unique.blocks.bed.tmp \
248.common_unique.blocks.bed.tmp \
249.common_unique.blocks.bed.tmp \
250.common_unique.blocks.bed.tmp \
251.common_unique.blocks.bed.tmp \
253.common_unique.blocks.bed.tmp \
254.common_unique.blocks.bed.tmp \
290.common_unique.blocks.bed.tmp \
291.common_unique.blocks.bed.tmp \
292.common_unique.blocks.bed.tmp \
|grep -v "all" >all.common_unique.blocks.bed
bedtools intersect -a all.common_unique.blocks.bed -b matloci.bed -nonamecheck> test
grep 'mito\|m86' all.common_unique.blocks.bed|cat test - >all.common_unique.blocks.matmito.bed
sed 's/PTT_//g' all.common_unique.blocks.bed|sed 's/0-1_contig_m86/mitochondria/g' |grep "chr\|mito" >all.common_unique.blocks.bed.tmp

mkdir Rplotfiles
cp all.common* Rplotfiles
cp ../all.* Rplotfiles
echo "################# Segregation Distortion Files #####################"
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

echo "#Remove all “failed” sequences --less than 50% coverage"
cut -f1-3 header|cat - ranges.bed| paste - \
2.ranges \
3.ranges \
4.ranges \
5.ranges \
6.ranges \
7.ranges \
9.ranges \
10.ranges \
11.ranges \
12.ranges \
13.ranges \
14.ranges \
15.ranges \
16.ranges \
17.ranges \
18.ranges \
20.ranges \
21.ranges \
22.ranges \
23.ranges \
24.ranges \
27.ranges \
28.ranges \
29.ranges \
32.ranges \
34.ranges \
35.ranges \
36.ranges \
37.ranges \
38.ranges \
41.ranges \
42.ranges \
43.ranges \
45.ranges \
46.ranges \
49.ranges \
50.ranges \
51.ranges \
52.ranges \
55.ranges \
65.ranges \
67.ranges \
68.ranges \
70.ranges \
71.ranges \
72.ranges \
73.ranges \
74.ranges \
75.ranges \
76.ranges \
77.ranges \
81.ranges \
85.ranges \
87.ranges \
88.ranges \
89.ranges \
90.ranges \
91.ranges \
92.ranges \
94.ranges \
95.ranges \
96.ranges \
99.ranges \
101.ranges \
102.ranges \
105.ranges \
108.ranges \
201.ranges \
203.ranges \
204.ranges \
205.ranges \
206.ranges \
207.ranges \
208.ranges \
210.ranges \
212.ranges \
214.ranges \
215.ranges \
216.ranges \
217.ranges \
219.ranges \
222.ranges \
226.ranges \
227.ranges \
229.ranges \
230.ranges \
233.ranges \
236.ranges \
238.ranges \
239.ranges \
240.ranges \
241.ranges \
242.ranges \
243.ranges \
244.ranges \
245.ranges \
246.ranges \
248.ranges \
249.ranges \
250.ranges \
251.ranges \
253.ranges \
254.ranges \
290.ranges \
291.ranges \
292.ranges \
>all.recomblocks.bed
sed 's/PTT_//g' all.recomblocks.bed |sed 's/0-1_contig_m86/mitochondria/g'  >all.recomblocks.bed.tmp

paste \
2.merge.ranges \
3.merge.ranges \
4.merge.ranges \
5.merge.ranges \
6.merge.ranges \
7.merge.ranges \
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
20.merge.ranges \
21.merge.ranges \
22.merge.ranges \
23.merge.ranges \
24.merge.ranges \
27.merge.ranges \
28.merge.ranges \
29.merge.ranges \
32.merge.ranges \
34.merge.ranges \
35.merge.ranges \
36.merge.ranges \
37.merge.ranges \
38.merge.ranges \
41.merge.ranges \
42.merge.ranges \
43.merge.ranges \
45.merge.ranges \
46.merge.ranges \
49.merge.ranges \
50.merge.ranges \
51.merge.ranges \
52.merge.ranges \
55.merge.ranges \
65.merge.ranges \
67.merge.ranges \
68.merge.ranges \
70.merge.ranges \
71.merge.ranges \
72.merge.ranges \
73.merge.ranges \
74.merge.ranges \
75.merge.ranges \
76.merge.ranges \
77.merge.ranges \
81.merge.ranges \
85.merge.ranges \
87.merge.ranges \
88.merge.ranges \
89.merge.ranges \
90.merge.ranges \
91.merge.ranges \
92.merge.ranges \
94.merge.ranges \
95.merge.ranges \
96.merge.ranges \
99.merge.ranges \
101.merge.ranges \
102.merge.ranges \
105.merge.ranges \
108.merge.ranges \
201.merge.ranges \
203.merge.ranges \
204.merge.ranges \
205.merge.ranges \
206.merge.ranges \
207.merge.ranges \
208.merge.ranges \
210.merge.ranges \
212.merge.ranges \
214.merge.ranges \
215.merge.ranges \
216.merge.ranges \
217.merge.ranges \
219.merge.ranges \
222.merge.ranges \
226.merge.ranges \
227.merge.ranges \
229.merge.ranges \
230.merge.ranges \
233.merge.ranges \
236.merge.ranges \
238.merge.ranges \
239.merge.ranges \
240.merge.ranges \
241.merge.ranges \
242.merge.ranges \
243.merge.ranges \
244.merge.ranges \
245.merge.ranges \
246.merge.ranges \
248.merge.ranges \
249.merge.ranges \
250.merge.ranges \
251.merge.ranges \
253.merge.ranges \
254.merge.ranges \
290.merge.ranges \
291.merge.ranges \
292.merge.ranges \
>all.recomblocks.vcf.tmp
rm *merge.ranges

echo "####################################################### Recombination Breakpoints ##################################################################"
cat \
2.break \
3.break \
4.break \
5.break \
6.break \
7.break \
9.break \
10.break \
11.break \
12.break \
13.break \
14.break \
15.break \
16.break \
17.break \
18.break \
20.break \
21.break \
22.break \
23.break \
24.break \
27.break \
28.break \
29.break \
32.break \
34.break \
35.break \
36.break \
37.break \
38.break \
41.break \
42.break \
43.break \
45.break \
46.break \
49.break \
50.break \
51.break \
52.break \
55.break \
65.break \
67.break \
68.break \
70.break \
71.break \
72.break \
73.break \
74.break \
75.break \
76.break \
77.break \
81.break \
85.break \
87.break \
88.break \
89.break \
90.break \
91.break \
92.break \
94.break \
95.break \
96.break \
99.break \
101.break \
102.break \
105.break \
108.break \
201.break \
203.break \
204.break \
205.break \
206.break \
207.break \
208.break \
210.break \
212.break \
214.break \
215.break \
216.break \
217.break \
219.break \
222.break \
226.break \
227.break \
229.break \
230.break \
233.break \
236.break \
238.break \
239.break \
240.break \
241.break \
242.break \
243.break \
244.break \
245.break \
246.break \
248.break \
249.break \
250.break \
251.break \
253.break \
254.break \
290.break \
291.break \
292.break \
|bedtools sort -i stdin -faidx $refPTT.fai \
>all.breaks

bedtools makewindows -g $refPTT.fai -w 10000 -s 10000 > genome.windows.bed
bedtools intersect -a genome.windows.bed -b all.breaks -c -nonamecheck \
|bedtools intersect -nonamecheck -wa -wb -a stdin -b all.breaks \
|awk '{print $0 "\t" $7-$6}'  > all.breaks.density

echo "####################################################### LD: r2 and r ##################################################################"
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

Create a r2 and r matrix 
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
Rscript --vanilla LinkageDisequilibrium.v2.R all.r.ptm.matrix.ld all.r2.ptm.matrix.ld all.r2.ptm.ld all.r.ptm.inter.ld
Rscript --vanilla LinkageDisequilibrium.v2.R all.r.ptt.matrix.ld all.r2.ptt.matrix.ld all.r2.ptt.ld all.r.ptt.inter.ld





