# Study of genetic incompatibilities in the plant pathogen Pyrenophora teres using lab crosses and field populations
scripts for analyzing Pyrenophora teres subspecies (PTT and PTM). 
## 1. 0-1 PTT and FGOBPtm-1 PTM Synteny of orthologous genes
### a. Circos Plot
Requires orthologous "all genes", "biosynthetic gene clusters", "effectors", "synteny chords" files
```
CirclePlot.bash.R
```
## 2. Lab Crosses
### a. Fitness of parents, intraspecies, and interspecies progeny
Fitness data is stored in PterespopsonKombar.csv
```
Pteres_cross_phenotypes_allrep.R
```
### b. Process raw sequencing data
Cleans and filters reads
Illumina and Ion Torrent sequences for PttxPtt, PtmxPtm, and PttxPtm progeny must be downloaded from NCBI-SRA
Reference genomes of 0-1, 15A, SG1, FGOBPtm-1 can also be found on NCBI
```
SL_PTTxPTT.parallel.v2.sh
SL_PTMxPTM.parallel.v2.sh
SL_PTTxPTMhybrid.parallel.v2.sh
```
Identifies recombination blocks in the PttxPtt, PtmxPtm, and PttxPtm progeny
Segregation Distortion files are also generated
```
PTTxPTTcross.sbatch.command.sh
PTMxPTMcross.sbatch.command.sh
PTTxPTMhybrid.sbatch.command.v2.sh
```
### c. Segregation Distortion
Visualizes Recombination Blocks, and Segregation Distortion with respect to genomic features (gene density, repeats, effectors, biosynthetic gene clusters, and Mat loci)
Calculates Recombination Hotspots
```
SegregationDistortion.R
```
### d. Linkage Disequilibirum
Calculates the odds ratio of parental vs hybrid genotypes from interchromosomal Linkage Disequilibirum 
```
LinkageDisequilibrium.v2.R
```
Example command reads in r and r2 as a matrix, and as a dataframe
```
Rscript --vanilla LinkageDisequilibrium.v2.R all.r.ptm.matrix.ld all.r2.ptm.matrix.ld all.r2.ptm.ld all.r.ptm.inter.ld
```
## 3. Field Population
### a. Cleans and filters reads
Illumina reads of Ptt and Ptm isolates from the Iranian population can be found on NCBI-SRA
```
SL_Populations_aligntoPTT.parallel.v2.sh
```
### b. Divergence Statistics (Dxy and Fst)
Calls SNPs in the the Ptt and Ptm Iranian population
```
SL_Populations_aligntoPTT.joint.v3.sh
```
Uses a sliding window to calculate Dxy and Fst
Requires scripts from https://github.com/simonhmartin/genomics_general
```
SL_DivergenceGeneFlow.v2.sh
```
Generates Figure for divergence statistics
```
SL_DivergenceGeneFlow.v2.sh
```
