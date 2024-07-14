## Extracting bfiles for filtered variants	[bfile_imp.sh]

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=240:0:0
#$ -l h_vmem=50G

module load plink/2.0-20200328

gdir="/data/Wolfson-PNU-dementia/imputed_genotypes"
ddir="/data/Wolfson-PNU-dementia/lungdev/bfile_imp_filt"

for i in {1..22}; do
plink2 \
 --bgen ${gdir}/ukb_imp_chr${i}_v3.bgen ref-first \
 --sample ${gdir}/ukb59138_imp_chr${i}_v3_s487296.sample \
 --extract variants_filtered_ID.txt \
 --make-bed \
 --out ${ddir}/ichr${i}_imp_filt
done