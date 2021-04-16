#!/bin/bash

#
#  submit by  sbatch Haplotype_caller.sh
#
#  specify the job name
#SBATCH --job-name=PTTxPTMhybrid
#  how many cpus are requested
#SBATCH --ntasks=4
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=90:00:00
#  maximum requested memory
#SBATCH --mem=70G
#  write std out and std error to these files
#SBATCH --error=PTTxPTMhybrid.%J.err
#SBATCH --output=PTTxPTMhybrid.%J.out
#  send a mail for job start, end, fail, etc.
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

########################################################################
refPTT=~/references/hybrid_references/PTT_0-1_assembly.v14.fa
refPTM=~/references/hybrid_references/PTM_FGOB10Ptm-1_assembly.v7.fasta
snps=~/PTTxPTMtest/SNPcalling/
parsnp_qry=~/PTTxPTMtest/SNPcalling/parsnp_qry/
parsnp_out=~/PTTxPTMtest/SNPcalling/parsnp_out/

#run bash to submit parallel jobs on the cluster
### example: sbatch SL_PTTxPTMhybrid.parallel.sh 100 /home/yuzon/PTTxPTMtest/
########################################################################
cd SNPcalling
cat \
1.PTMPTT.flt.bed \
2.PTMPTT.flt.bed \
3.PTMPTT.flt.bed \
4.PTMPTT.flt.bed \
5.PTMPTT.flt.bed \
6.PTMPTT.flt.bed \
7.PTMPTT.flt.bed \
8.PTMPTT.flt.bed \
9.PTMPTT.flt.bed \
10.PTMPTT.flt.bed \
11.PTMPTT.flt.bed \
12.PTMPTT.flt.bed \
13.PTMPTT.flt.bed \
14.PTMPTT.flt.bed \
15.PTMPTT.flt.bed \
16.PTMPTT.flt.bed \
17.PTMPTT.flt.bed \
18.PTMPTT.flt.bed \
19.PTMPTT.flt.bed \
20.PTMPTT.flt.bed \
21.PTMPTT.flt.bed \
22.PTMPTT.flt.bed \
23.PTMPTT.flt.bed \
24.PTMPTT.flt.bed \
25.PTMPTT.flt.bed \
26.PTMPTT.flt.bed \
27.PTMPTT.flt.bed \
28.PTMPTT.flt.bed \
29.PTMPTT.flt.bed \
30.PTMPTT.flt.bed \
31.PTMPTT.flt.bed \
32.PTMPTT.flt.bed \
33.PTMPTT.flt.bed \
34.PTMPTT.flt.bed \
35.PTMPTT.flt.bed \
36.PTMPTT.flt.bed \
37.PTMPTT.flt.bed \
38.PTMPTT.flt.bed \
39.PTMPTT.flt.bed \
40.PTMPTT.flt.bed \
41.PTMPTT.flt.bed \
42.PTMPTT.flt.bed \
43.PTMPTT.flt.bed \
44.PTMPTT.flt.bed \
45.PTMPTT.flt.bed \
46.PTMPTT.flt.bed \
47.PTMPTT.flt.bed \
48.PTMPTT.flt.bed \
49.PTMPTT.flt.bed \
50.PTMPTT.flt.bed \
51.PTMPTT.flt.bed \
52.PTMPTT.flt.bed \
53.PTMPTT.flt.bed \
54.PTMPTT.flt.bed \
55.PTMPTT.flt.bed \
56.PTMPTT.flt.bed \
57.PTMPTT.flt.bed \
58.PTMPTT.flt.bed \
59.PTMPTT.flt.bed \
60.PTMPTT.flt.bed \
61.PTMPTT.flt.bed \
62.PTMPTT.flt.bed \
63.PTMPTT.flt.bed \
64.PTMPTT.flt.bed \
65.PTMPTT.flt.bed \
66.PTMPTT.flt.bed \
67.PTMPTT.flt.bed \
68.PTMPTT.flt.bed \
69.PTMPTT.flt.bed \
70.PTMPTT.flt.bed \
72.PTMPTT.flt.bed \
73.PTMPTT.flt.bed \
74.PTMPTT.flt.bed \
75.PTMPTT.flt.bed \
76.PTMPTT.flt.bed \
77.PTMPTT.flt.bed \
78.PTMPTT.flt.bed \
79.PTMPTT.flt.bed \
79b.PTMPTT.flt.bed \
80.PTMPTT.flt.bed \
81.PTMPTT.flt.bed \
82.PTMPTT.flt.bed \
83.PTMPTT.flt.bed \
84.PTMPTT.flt.bed \
85.PTMPTT.flt.bed \
86.PTMPTT.flt.bed \
87.PTMPTT.flt.bed \
88.PTMPTT.flt.bed \
89.PTMPTT.flt.bed \
90.PTMPTT.flt.bed \
91.PTMPTT.flt.bed \
92.PTMPTT.flt.bed \
93.PTMPTT.flt.bed \
94.PTMPTT.flt.bed \
95.PTMPTT.flt.bed \
96.PTMPTT.flt.bed \
97.PTMPTT.flt.bed \
98.PTMPTT.flt.bed \
99.PTMPTT.flt.bed \
100.PTMPTT.flt.bed \
101.PTMPTT.flt.bed \
102.PTMPTT.flt.bed \
103.PTMPTT.flt.bed \
104.PTMPTT.flt.bed \
105.PTMPTT.flt.bed \
106.PTMPTT.flt.bed \
107.PTMPTT.flt.bed \
108.PTMPTT.flt.bed \
109.PTMPTT.flt.bed \
110.PTMPTT.flt.bed \
111.PTMPTT.flt.bed \
112.PTMPTT.flt.bed \
113.PTMPTT.flt.bed \
114.PTMPTT.flt.bed \
115.PTMPTT.flt.bed \
116.PTMPTT.flt.bed \
117.PTMPTT.flt.bed \
118.PTMPTT.flt.bed \
119.PTMPTT.flt.bed \
120.PTMPTT.flt.bed \
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


cat *.common_unique.blocks.bed.tmp |grep -v "all" >all.common_unique.blocks.bed
bedtools intersect -a all.common_unique.blocks.bed -b matloci.bed -nonamecheck> test
grep 'mito\|m86' all.common_unique.blocks.bed|cat test - >all.common_unique.blocks.matmito.bed

mkdir Rplotfiles
cp all.common* Rplotfiles
cp ../all.* Rplotfiles
################# Segregation Distortion Files #####################
##awk '$6==1 & $7==1 {print $0}' SNPcalling/nucmer/all.recomblocks.bed ##check if error such that PTT and PTM alles in an isolate overlap (begin with even=PTM alleles odd=PTT alleles)
cat *support|cut -f1,2| sort -k1,1 -k2,2n  >start
cat *support|cut -f1,3| sort -k1,1 -k2,2n  >stop
cat start stop | sort -k1,1 -k2,2n |uniq>tmp.ranges
tail -n +2 tmp.ranges|cut -f1,2|paste tmp.ranges -|awk '$1==$3 {print $1 "\t" $2 "\t" $4}'|bedtools sort -i stdin -faidx $refPTT.fai >ranges.bed

for each in *PTM.support
do
entry=${each%.PTM.support}
bedtools intersect -a ranges.bed -b $entry.PTM.support -wao -nonamecheck > $entry.PTM.ranges.bed
bedtools intersect -a ranges.bed -b $entry.PTT.support -wao -nonamecheck > $entry.PTT.ranges.bed
echo -e "chrom\tstart\t stop\t$entry.PTM.support\t$entry.PTT.support" >header
paste  $entry.PTM.ranges.bed  $entry.PTT.ranges.bed |cut -f1-3,7,17 | cat header - |cut -f4,5> $entry.ranges
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
1.ranges \
2.ranges \
3.ranges \
5.ranges \
6.ranges \
7.ranges \
8.ranges \
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
19.ranges \
20.ranges \
21.ranges \
22.ranges \
23.ranges \
25.ranges \
27.ranges \
28.ranges \
29.ranges \
30.ranges \
31.ranges \
32.ranges \
33.ranges \
34.ranges \
36.ranges \
37.ranges \
38.ranges \
39.ranges \
40.ranges \
41.ranges \
42.ranges \
43.ranges \
45.ranges \
46.ranges \
47.ranges \
48.ranges \
49.ranges \
50.ranges \
51.ranges \
52.ranges \
53.ranges \
55.ranges \
56.ranges \
57.ranges \
58.ranges \
59.ranges \
60.ranges \
61.ranges \
62.ranges \
63.ranges \
64.ranges \
65.ranges \
66.ranges \
67.ranges \
68.ranges \
69.ranges \
70.ranges \
72.ranges \
75.ranges \
76.ranges \
77.ranges \
79b.ranges \
81.ranges \
82.ranges \
83.ranges \
84.ranges \
85.ranges \
86.ranges \
87.ranges \
88.ranges \
90.ranges \
91.ranges \
92.ranges \
93.ranges \
95.ranges \
96.ranges \
97.ranges \
98.ranges \
100.ranges \
101.ranges \
102.ranges \
103.ranges \
104.ranges \
105.ranges \
106.ranges \
107.ranges \
108.ranges \
109.ranges \
110.ranges \
111.ranges \
112.ranges \
113.ranges \
114.ranges \
115.ranges \
116.ranges \
117.ranges \
118.ranges \
119.ranges \
120.ranges \
>all.recomblocks.bed

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

####GENETIC MAP: mapthis dataframe####
cat  all.recomblocks.vcf.tmp |tail -n +2|sed 's/\./-/g'|sed 's/PTM/B/g' |sed 's/PTT/A/g' >matrix
head -n 1 all.recomblocks.vcf.tmp  >matrix_header
cut -f1,2,3 all.recomblocks.bed |awk -v OFS="\t" '{print $1 ":" $2"-"$3, $1}'| awk -v OFS="\t" '{gsub(/chrPTT_/,"",$2); print}'>matrix_col
cat matrix_header matrix|paste matrix_col - | datamash transpose >test

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
1.break \
2.break \
3.break \
5.break \
6.break \
7.break \
8.break \
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
19.break \
20.break \
21.break \
22.break \
23.break \
25.break \
27.break \
28.break \
29.break \
30.break \
31.break \
32.break \
33.break \
34.break \
36.break \
37.break \
38.break \
39.break \
40.break \
41.break \
42.break \
43.break \
45.break \
46.break \
47.break \
48.break \
49.break \
50.break \
51.break \
52.break \
53.break \
55.break \
56.break \
57.break \
58.break \
59.break \
60.break \
61.break \
62.break \
63.break \
64.break \
65.break \
66.break \
67.break \
68.break \
69.break \
70.break \
72.break \
75.break \
76.break \
77.break \
79b.break \
81.break \
82.break \
83.break \
84.break \
85.break \
86.break \
87.break \
88.break \
90.break \
91.break \
92.break \
93.break \
95.break \
96.break \
97.break \
98.break \
100.break \
101.break \
102.break \
103.break \
104.break \
105.break \
106.break \
107.break \
108.break \
109.break \
110.break \
111.break \
112.break \
113.break \
114.break \
115.break \
116.break \
117.break \
118.break \
119.break \
120.break \
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




