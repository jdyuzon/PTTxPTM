# Study of genetic incompatibilities in the plant pathogen Pyrenophora teres using lab crosses and field populations 
scripts for analyzing Pyrenophora teres subspecies (PTT and PTM). 
## 1. Lab Crosses
### a. Fitness of parents, intraspecies, and interspecies progeny
```
Pteres_cross_phenotypes_allrep.R
```
### b. Process raw sequencing data
Cleans and filters reads
```
SL_PTTxPTT.parallel.v2.sh
SL_PTMxPTM.parallel.v2.sh
SL_PTTxPTMhybrid.parallel.v2.sh
```
Identifies recombination blocks
```
PTTxPTTcross.sbatch.command.sh
PTMxPTMcross.sbatch.command.sh
PTTxPTMhybrid.sbatch.command.v2.sh
```
### c. Segregation Distortion
```
SegregationDistortion.R
```
### d. Linkage Disequilibirum
```
LinkageDisequilibrium.v2.R
```
## 2. Plant Pathogen Genome Compartmentalization
```
SegregationDistortion_twospeed.R
```
## 3. Field Population
### a. Cleans and filters reads
```
SL_Populations_aligntoPTT.parallel.v2.sh
```
### b. Divergence Statistics (Dxy and Fist)
requires scripts from https://github.com/simonhmartin/genomics_general
```
SL_Populations_aligntoPTT.joint.v3.sh
SL_DivergenceGeneFlow.v2.sh
```
