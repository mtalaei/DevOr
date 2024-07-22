#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=240:0:0
#$ -l h_vmem=40G
#$ -l highmem

module load R
Rscript Extract_CogDem.R

