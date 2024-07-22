#-------------------------------------------------------------------------------
# Lung Function data management (instance 0 @Baseline) - As an outcome for Regenie
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) 

library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")
library(naniar)
library(ukbtools)


##====Read  (^)
#raw extracted file
DB <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_Cov_LF.rds")
str(DB)
#managed file
DBs <- readRDS("/data/Wolfson-PNU-dementia/lungdev/UKB_Projects/ukb_pheno_LF.rds")


##====Selected variables at instances 0 
DBs <- subset(DB, select=c(FID, IID, FEV1_b0, FVC_b0, FEV1FVC_b0))
names(DBs)

#=======================================================================Exclusions
#====Missing structure: Original
map(DBs, ~sum(!is.na(.)))

#====Ethnicity from CoreDataClin1 to exclude non-Whites
DB_Cov2 <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_ethn.rds")
table(DB_Cov2$Ethn_b)                              #0 472658, 1 28905
map(DB_Cov2, ~sum(is.na(.)))
#merge ethnicity 
dim(DBs)                                        #502462     5
DBs <- left_join(DBs,DB_Cov2,by="IID")
names(DBs)
dim(DBs)                                        #502462     7
DBs %>% count(Ethn_b)                           #0 472658, 1 28905, NA 899
#exclude non-Whites
DBs=subset(DBs, Ethn_b==0 & !is.na(Ethn_b)) 
DBs %>% count(Ethn_b)                           #0 472658
DBs=subset(DBs, select= -c(Ethn,Ethn_b))
names(DBs)

#Missing structure: non-Whites
map(DBs, ~sum(!is.na(.)))


#=====Kinship (Relatedness)
kinship = read_table("/data/Wolfson-PNU-dementia/UKB/ukb_helper_files/ukb59138_rel_s488264.dat")
summary(kinship$Kinship)
#Min.    1st Qu. Median  Mean    3rd Qu. Max. 
#0.0442  0.0566  0.0699  0.1185  0.2254  0.4999         n=? 107116
kinship$Kinship_0.125 <- ifelse(kinship$Kinship >= 0.125, c("higher"), c("lower"))
table(kinship$Kinship_0.125)                            #n=73090/34026
#1st degree relatives (kinship coefficient>=0.125)
Exclu_kinship = kinship %>% filter(Kinship>0.125) %>% select(ID1) %>% rename("IID"="ID1") 
#Exclusion/Filtering
DBs <- DBs  %>% mutate(exclu_Kinship = ifelse(IID %in% Exclu_kinship$IID, "1", "0"))
table(DBs$exclu_Kinship)                              #n=442606  30052
#DBs<-filter(DBs, exclu_Kinship==0)
DBs=subset(DBs, exclu_Kinship==0)
table(DBs$exclu_Kinship)                              #n=442606
DBs=subset(DBs, select= -c(exclu_Kinship))
names(DBs)

#Missing structure: non-1st relatives
map(DBs, ~sum(!is.na(.)))


#=====Poor heterozygosity/missingness
#Read ukb_SQC_outliers.rds
saveRDS(DB_hmid, file="/data/Wolfson-PNU-dementia/lungdev/ukb_SQC_outliers.rds")
DB_hmid <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_SQC_outliers.rds")
table(DB_hmid$het.missing.outliers)                #968
DB_hmid=subset(DB_hmid, select= -c(V25))

#Merge het.missing.outliers to the outcome dataset
dim(DBs)                                        #442563     5
#.DBs <- merge(DBs,DB_hmid,by="IID")
DBs <- left_join(DBs,DB_hmid,by="IID")
names(DBs)
dim(DBs)                                        #442563     6
DBs %>% count(het.missing.outliers)                           #0 429138, 1 907, NA 12561
#exclude outliers
DBs=subset(DBs, het.missing.outliers==0 | is.na(het.missing.outliers)) 
DBs %>% count(het.missing.outliers)                           #0 429138, NA 12561
DBs=subset(DBs, select= -c(het.missing.outliers))
names(DBs)

#Missing structure: No poor heterozygosity/missingness
map(DBs, ~sum(!is.na(.)))


#=====Covariates: PC 
##Read Cov_LF file
DBc <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")
DBcs=subset(DBc, select= c(IID,age,sex,centre,pc1,pc5, Smoking_0, Height_0))
map(DBcs, ~sum(!is.na(.)))

DBcs <- DBcs %>% mutate(PCmis = ifelse(!is.na(pc1) & !is.na(pc5) & !is.na(age),0,1))
table(DBcs$PCmis)                               #488221  14241
DBcss=subset(DBcs, select= c(IID,PCmis))
str(DBcss)

#Merge Covariates&LF missing to the outcome dataset
dim(DBs)                                        #441699     5
DBs <- left_join(DBs,DBcss,by="IID")
names(DBs)
dim(DBs)                                        #441699     6
#exclude those without genetic data
DBs %>% count(PCmis)                            #0 429138, 1 12561
DBs=subset(DBs, PCmis==0) 
DBs %>% count(PCmis)                            #0 429138
DBs=subset(DBs, select= -c(PCmis))
names(DBs)

#Missing structure: No PC missing
map(DBs, ~sum(!is.na(.)))


#=====other Covariates: Smoking & Height
#Smoking
table(DBcs$Smoking_0)
DBcs <- DBcs %>% mutate(Smmis = ifelse(Smoking_0==-3,0,1))
table(DBcs$Smmis)                               #499513  2057
#Height
DBcs <- DBcs %>% mutate(Htmis = ifelse(!is.na(Height_0),0,1))
table(DBcs$Htmis)                               #499923  2539

DBcss=subset(DBcs, select= c(IID,Smmis,Htmis))
str(DBcss)

#Merge Covariates missing to the outcome dataset
dim(DBs)                                        #429138     5
DBs <- left_join(DBs,DBcss,by="IID")
names(DBs)
dim(DBs)                                        #429138     7
#exclude those without smoking data
DBs %>% count(Smmis)                            #0 1537, 1 427601
DBs=subset(DBs, Smmis==1) 
DBs %>% count(Smmis)                            #1 427601
DBs=subset(DBs, select= -c(Smmis))
names(DBs)

#exclude FVC without height data
DBs %>% count(Htmis)                            #0 426644, 1 957
DBs=subset(DBs, Htmis==1) 
DBs$FVC_b0[DBs$Htmis==1] <- NA
DBs=subset(DBs, select= -c(Htmis))
names(DBs)

#Missing structure: No PC missing
map(DBs, ~sum(!is.na(.)))

#=====Saved as: ukb_pheno_cogTS.rds
print("Saved as: ukb_pheno_LF.rds")
saveRDS(DBs, file="/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_LF.rds")
print("Saved as: ukb_pheno_LF.txt")
write.table(DBs, file = "/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_LF.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





###============================================Outcome files
##Read extracted file, managed
DBs <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_LF.rds")

#1.FVC
FVC_0 = DBs%>% select(contains("FID"), contains("IID"), FVC_b0)
map(FVC_0, ~sum(!is.na(.)))
FVC_0 <- FVC_0[complete.cases(FVC_0),]
print("Structure- FVC best measure at instance 0:")
str(FVC_0)
write.table(FVC_0, file = "/data/Wolfson-PNU-dementia/lungdev/FVC_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: FVC_0.txt")

#1.Ratio
Ratio_0 = DBs%>% select(contains("FID"), contains("IID"), FEV1FVC_b0)
map(Ratio_0, ~sum(!is.na(.)))
Ratio_0 <- Ratio_0[complete.cases(Ratio_0),]
print("Structure- Ratio best measure at instance 0:")
str(Ratio_0)
write.table(Ratio_0, file = "/data/Wolfson-PNU-dementia/lungdev/Ratio_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: Ratio_0.txt")
