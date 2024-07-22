#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=30G
module load plink
ddir="/data/Wolfson-PNU-dementia/lungdev"
sdir="/data/Wolfson-PNU-dementia/genotype_hardcalls"

rm -f ${ddir}/list_beds.txt
for chr in {2..22}; do echo "${sdir}/ukb_cal_chr${chr}_v2.bed ${sdir}/ukb_snp_chr${chr}_v2.bim ${sdir}/ukb59138_cal_chr1_v2_s488264.fam" >> ${ddir}/list_beds.txt; done

plink \
  --bed ${sdir}/ukb_cal_chr1_v2.bed \
  --bim ${sdir}/ukb_snp_chr1_v2.bim \
  --fam ${sdir}/ukb59138_cal_chr1_v2_s488264.fam \
  --merge-list ${ddir}/list_beds.txt \
  --make-bed --out ${ddir}/ukb_cal_allChrs

