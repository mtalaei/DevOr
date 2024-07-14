## Creating LD matrix for SuSiE	[ld_susie.sh]

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=15G

module load plink/1.9-170906
ddir="/data/Wolfson-PNU-dementia/lungdev"

plink \
  --bfile ${ddir}/ukb_imp_allChrs_filt \
  --extract variants_filtered_ID.txt \
  --r square \
  --keep-allele-order \
  --write-snplist \
  --out ${ddir}/ld_matrix
echo "Number of columns"
awk '{print NF}' ${ddir}/ld_matrix.ld | sort -nu | tail -n 1
echo "Number of rows"
cat ${ddir}/ld_matrix.ld | wc -l