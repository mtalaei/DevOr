##Regenie-1  Ratio	(FEV1/FVC ratio at baseline)

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=240:0:0
#$ -l h_vmem=10G

##Note: Lines 13, 21, and 30 need to be modified for each outcome

ddir="/data/Wolfson-PNU-dementia/lungdev"         #Data directory
odir="/data/home/hmy431/outputs/Ratio_0"          #Output directory

./regenie_v2.2.4.gz_x86_64_Centos7_mkl  \
  --step 1 \
  --bed ${ddir}/ukb_cal_allChrs \
  --extract ${ddir}/qc_pass.snplist \
  --keep ${ddir}/qc_pass.id \
  --phenoFile ${ddir}/Ratio_0.txt \
  --phenoColList FEV1FVC_b0 \
  --covarFile ${ddir}/ukb_Cov_LF.txt \
  --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 \
  --catCovarList centre,array,Smoking_0 \
  --maxCatLevels 22 \
  --threads 8 \
  --bsize 1000 \
  --strict \
  --lowmem \
  --out ${odir}/ukb_step1_Ratio_0
