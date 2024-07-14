#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=1:0:0
#$ -l h_vmem=8G

module load R
Rscript Prep_Coloc_NM_0.R
