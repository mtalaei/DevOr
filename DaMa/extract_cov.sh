#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=20G
#$ -l highmem

module load R
Rscript Extract_CoV.R

