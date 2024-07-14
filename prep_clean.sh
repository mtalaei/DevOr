#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=15G
module load plink/2.0-20200328
ddir="/data/Wolfson-PNU-dementia/lungdev"

plink2 \
  --bfile ${ddir}/ukb_cal_allChrs \
  --maf 0.05 --mac 10 --geno 0.015 --hwe 1e-6 \
  --write-snplist --write-samples --no-id-header \
  --out ${ddir}/qc_pass

