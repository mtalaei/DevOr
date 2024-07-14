##Preparing Regenie output  *Reaction time	(Mean time to correctly identify matches.2.0) 
#For colocalization analysis (adding MAF, creating P-value, LF data, etc.)
#Note: Lines 11, 62, 64, and 66 need to be modified for each outcome

#----@HPC
library(dplyr)
library(tidyverse)
setwd("/data/home/hmy431/outputs/")     

#Adding MAF
DB1 <- read.table("gF9i_2/ukb_step2_gF9i_2.regenie", header=TRUE)
print("Dimension of output dataset:")
dim(DB1)
print("Duplicate in output dataset:")
sum(duplicated(DB1$ID))
DB2 <- read.table("/data/Wolfson-PNU-dementia/lungdev/ukb_maf_v3.txt", header=TRUE)
print("Dimension of UKB-mfi dataset:")
dim(DB2)
print("Duplicate in UKB-mfi dataset:")
sum(duplicated(DB2$ID))
DB2 <- rename(DB2, ALLELE0=Allele2)
DB2 <- rename(DB2, ALLELE1=Allele1)
DB <- merge(DB1, DB2, by=c("ID","ALLELE0","ALLELE1"))
print("Dimension of merged dataset:")
dim(DB)
print("Duplicate in merged dataset:")
sum(duplicated(DB$ID))
print("Structure of merged dataset:")
str(DB)

#Removing EXTRA 
print("EXTRA content")
table(DB$EXTRA)
summary(DB$EXTRA)
DB=subset(DB, select= -c(EXTRA))

#Calculating P-value
DB$LOG10P <- as.numeric(DB$LOG10P)
DB$pval <- 10^(DB$LOG10P*-1)
print("P-value calculated")

#Selecting variables needed for co-localization analysis
DB=subset(DB, select= c(ID, INFO, MAF, BETA, SE, pval))

#Adding further items
DB3 <- read.table("fromLaura/LF_statistics.txt", header=TRUE)
print("Dimension of LF dataset:")
dim(DB3)
DB4 <- read.table("fromLaura/variants_filtered_ID_gene.txt", header=TRUE)
print("Dimension of variants dataset:")
dim(DB4)
DB <- merge(DB, DB3, by="ID")
print("Dimension of merged-LF dataset:")
dim(DB4)
DB <- merge(DB, DB4, by="ID")
print("Dimension of merged-variants dataset:")
dim(DB)
print("Structure of final merged dataset:")
str(DB)

#Save
saveRDS(DB, file='gF9i_2/ukb_step2_gF9i_2_add.rds')
#Tab delimited text file
write.table(DB, file = "gF9i_2/ukb_step2_gF9i_2_add.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#Compressed file
write.table(DB, gzfile("gF9i_2/ukb_step2_gF9i_2_add.txt.gz"), row.names = FALSE, col.names = TRUE, quote = FALSE) 