##Regenie-1  Pairs matching adjusted for height	(Number of incorrect matches in round.0.1 -> R_Test1_PairsMatching_f399_0)

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=240:0:0
#$ -l h_vmem=10G

##Note: Lines 13, 21, and 30 need to be modified for each outcome

ddir="/data/Wolfson-PNU-dementia/lungdev"            #Data directory
odir="/data/home/hmy431/outputs/PairsMatching_0_H"    #Output directory

./regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 1 \
  --bed ${ddir}/ukb_cal_allChrs \
  --extract ${ddir}/qc_pass.snplist \
  --keep ${ddir}/qc_pass.id \
  --phenoFile ${ddir}/PairsMatching_0.txt \
  --phenoColList PairsMatching_0_logrz \
  --covarFile ${ddir}/ukb_Cov_LF.txt \
  --covarColList age,sex,Height_0,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 \
  --catCovarList centre,array \
  --maxCatLevels 22 \
  --threads 8 \
  --bsize 1000 \
  --strict \
  --lowmem \
  --out ${odir}/ukb_step1_Pairs_0
