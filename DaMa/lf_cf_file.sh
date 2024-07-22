#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=40G
#$ -l highmem

module load R
Rscript LF_CF_file.R

