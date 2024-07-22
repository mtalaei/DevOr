##Preparing file for colocalisation analysis using GWAS data
#For colocalization analysis (adding MAF, creating P-value, LF data, etc.)
#===========================================================
rm(list = ls(all.names = TRUE)) 


#----@HPC
library(dplyr)
library(tidyverse)
library(data.table)
library(R.utils)
library("tidyr")
setwd("/data/home/hmy431/")     


##Summary statistics for Alzheimer’s disease from Bellenguez et al., 2022
#https://www.ebi.ac.uk/gwas/publications/35379992
#http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz 

#====Reading downloaded GWAS summary data
df <- fread('/data/home/hmy431/data/GCST90027158_buildGRCh38.tsv.gz')
str(df)                                           #21,101,114 obs. of  17 variables
df$pval <- as.numeric(df$p_value)
dfs <- subset(df, pval<5e-8)                      #5,637 SNPs with P<5*10-8
df <- rename(df, ID=variant_id)
names(df)

#split variant_alternate_id
df <- df %>% separate(variant_alternate_id, c("chromosome","position", "reference_allele", "alternate_allele"), ":")    
str(df)                                           #21,101,114 obs. of  20 variables
sum(duplicated(df$ID))                            #30,648 Duplicates in GWAS dataset
sum(duplicated(df$variant_alternate_id))          #0
df <- rename(df, Allele1=reference_allele)
df <- rename(df, Allele2=alternate_allele)

#==Adding MAF
DB_MAF <- read.table("/data/Wolfson-PNU-dementia/lungdev/ukb_maf_v3.txt", header=TRUE)
str(DB_MAF)                                           #93,095,623 obs. of  4 variables
sum(duplicated(DB_MAF$ID))      #320,321 Duplicates in UKB-mfi dataset
#Merge based on ID
#DB <- merge(df, DB_MAF, by="ID")                    
#dim(DB)                                             #19,265,975
#Merge based on ID, Allele1, and Allele2
DB <- merge(df, DB_MAF, by=c("ID","Allele1","Allele2"))
str(DB)                                             #19,067,815
#sum(duplicated(DB_MAF$ID, DB_MAF$Allele1)) #?


#==Selecting variables
#DBs=subset(DB, select= c(ID, INFO, MAF, BETA, SE, pval))   #Selecting variables needed for co-localization analysis
DBs=subset(DB, select= c(ID, effect_allele, other_allele, Allele1, Allele2, beta, 
                         standard_error, pval, MAF, n_cases, n_controls, base_pair_location))
DBs <- rename(DBs, BETA=beta, SE=standard_error, POS=base_pair_location)
str(DBs)

#====Adding further items
#==Variants dataset
DB_V <- read.table("outputs/fromLaura/variants_filtered_ID_gene.txt", header=TRUE)
str(DB_V)                                           #15,640
#Merge with variants dataset
DBs <- merge(DBs, DB_V, by="ID")                    
str(DBs)                                             #12,888 

#==LF dataset:
DB_LF <- read.table("outputs/fromLaura/LF_statistics.txt", header=TRUE)
str(DB_LF)                                          #15,298
DB_LFs <- subset(DB_LF, p_fvc<5e-8 | p_ratio<5e-8)
dim(DB_LFs)                                         #709 SNPs with P<5*10-8
#Merge with LF dataset
DBs <- merge(DBs, DB_LF, by="ID")
str(DBs)                                             #12,598 

sum(duplicated(DBs$ID))                              #1 Duplicate in merged dataset
DBs <- DBs %>% distinct(ID, .keep_all= TRUE)
str(DBs)                                             #12,597 


#Prepare original extracted variants (no selection or addition)
DBso <- DB  %>% mutate(DevOr = ifelse(ID %in% DB_V$ID, "1", "0"))
table(DBso$DevOr)                    #0  19054927      1 12888
DBso=subset(DBso, DevOr==1)
str(DBso)


#====Save
#original extraction: includes all variables available in the source 
saveRDS(DBso, file='outputs/GWAS_AD/ukb_GWASAD.rds')

#final version for analysis
saveRDS(DBs, file='outputs/GWAS_AD/ukb_GWASAD_add.rds')
#Tab delimited text file
write.table(DBs, file = "outputs/GWAS_AD/ukb_GWASAD_add.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#Compressed file
write.table(DBs, gzfile("outputs/GWAS_AD/ukb_GWASAD_add.txt.gz"), row.names = FALSE, col.names = TRUE, quote = FALSE) 
