## Merging extracted bfiles for filtered variants	[bfile_merge.sh]

#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:0:0
#$ -l h_vmem=10G

module load plink/1.9-170906


ddir="/data/Wolfson-PNU-dementia/lungdev"

plink \
 --bfile ${ddir}/bfile_imp_filt/ichr1_imp_filt_tmp \
 --merge-list ${ddir}/bfile_imp_filt/allfiles_tmp.txt \
 --make-bed \
 --out ${ddir}/ukb_imp_allChrs_filt
