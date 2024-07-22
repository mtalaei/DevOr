#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=1:0:0
#$ -l h_vmem=8G

module load R
echo "---------------------Fluid intelligence/reasoning (Fluid intelligence score.2.0)"
echo "---------------------Prep_Coloc_FI.R"
Rscript Prep_Coloc_FI_2.R
