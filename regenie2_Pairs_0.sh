##Regenie-2  Pairs Matching	(Number of incorrect matches in round.0.1 -> R_Test1_PairsMatching_f399_0)

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=240:0:0
#$ -l h_vmem=10G

##Note: Lines 15, 24, 30, 39, 41, and 43 need to be modified for each outcome

##Defining directories
ddir="/data/Wolfson-PNU-dementia/lungdev"             #Data directory
gdir="/data/Wolfson-PNU-dementia/imputed_genotypes"   #Genetic data directory
odir="/data/home/hmy431/outputs/PairsMatching_0"     #Output directory

##Running Regenie-step2 for each chromosome
for i in {1..22}; do
./regenie_v2.2.4.gz_x86_64_Centos7_mkl  \
  --step 2 \
  --bgen  ${gdir}/ukb_imp_chr${i}_v3.bgen \
  --sample ${gdir}/ukb59138_imp_chr${i}_v3_s487296.sample \
  --phenoFile ${ddir}/PairsMatching_0.txt \
  --phenoColList PairsMatching_0_logrz \
  --covarFile ${ddir}/ukb_cov.txt \
  --covarColList age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 \
  --catCovarList centre,array \
  --maxCatLevels 22 \
  --extract variants_filtered_ID.txt \
  --pred ${odir}/ukb_step1_Pairs_0_pred.list \
  --bsize 1000 \
  --threads 8 \
  --strict \
  --minMAC 10 \
  --out ${odir}/ukb_step2_chr$i
done

##Appending output files (without a repeated first row)
cp ${odir}/ukb_step2_chr1_PairsMatching_0_logrz.regenie ${odir}/ukb_step2_Pairs_0.regenie
for i in {2..22}; do
cp ${odir}/ukb_step2_chr${i}_PairsMatching_0_logrz.regenie ${odir}/ukb_step2_chr${i}
sed -i '1d' ${odir}/ukb_step2_chr${i}
cat ${odir}/ukb_step2_chr${i} >> ${odir}/ukb_step2_Pairs_0.regenie
rm ${odir}/ukb_step2_chr${i}
done