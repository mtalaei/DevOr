#-------------------------------------------------------------------------------
# Dementia data management 
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) 

library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")
library(naniar)
#library(plyr)

##====Read 
#raw extracted file
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/UKBW_52532_110822A_CogDem.rds")

#managed file
DBs <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_Dem.rds")

#old
#DB <- read_csv(file="/data/Wolfson-PNU-dementia/lungdev/gfactor_and_updateddiag/NewCinicalDefs.csv")
#DBs = DB%>% select(contains("EID"), contains("Dementia"), contains("Alzheimers"))
#saveRDS(DBs, file="/data/Wolfson-PNU-dementia/lungdev/NewCinicalDefs_Dementia.rds")


#=====Outcome selection
#Note: pDementia_0 (ML_C42C240Xf41270_DementiaDiagyrBin0) has later replaced with a more accurate variable due to improvement in time variable
DBs=subset(DB, select= c("EID", 
                         "ML_C42C240Xf41270_Alzheimers", 
                         "ML_C42C240Xf41270_VascularDementia", 
                         "ML_C42C240Xf41270_Dementia", 
                         "ML_C42C240Xf41270_DementiaTenure",
                         "ML_C42C240Xf41270_DementiaDiagyrBin0",
                         "ML_C42C240Xf41270_DementiaAge",
                         "R_Age_f21003_0",
                         "R_LTFU_f191")) 
print("Structure of selected Dementia variables:")
str(DBs)

DBs <- rename(DBs, c(IID=EID))
DBs$FID <- DBs$IID
DBs <- DBs %>% relocate(FID, .before = IID)
DBs <- rename(DBs,
              DemAlz=ML_C42C240Xf41270_Alzheimers,
              DemVas=ML_C42C240Xf41270_VascularDementia,
              Dementia=ML_C42C240Xf41270_Dementia,
              DementiaTenure=ML_C42C240Xf41270_DementiaTenure,
              pDementia_0=ML_C42C240Xf41270_DementiaDiagyrBin0,
              DementiaAge=ML_C42C240Xf41270_DementiaAge,
              Age_0=R_Age_f21003_0,
              LTFU_0=R_LTFU_f191)
str(DBs)              


describe(DBs$DementiaAge)           #min 38.9 max 84.41 n 7382
describe(DBs$Age_0)                 #min 37 max 73  n 502411   
describe(DBs$DementiaTenure)        #min -12.71 max 16.19  n 501771   
#DBs <- DBs %>% mutate(pDementia_02 = ifelse(DementiaAge<=Age_0,1,0))
table(DBs$pDementia_0)              #502299    112
#table(DBs$pDementia_02)             #7270  112

#Dementia frequency
print("Frequency of Alzheimer Dementia:")
table(DBs$DemAlz)
print("Frequency of Vascular Dementia:")
table(DBs$DemVas)
print("Frequency of all Dementia:")
table(DBs$Dementia)


##Set working directory for outputs (figures)
setwd("/data/home/hmy431/outputs/Figures/")
#Missing structure
map(DBs, ~sum(!is.na(.)))


#====Ethnicity from CoreDataClin1 to exclude non-Whites
DB_Cov2 <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_ethn.rds")
table(DB_Cov2$Ethn_b)                              #0 472658, 1 28905
map(DB_Cov2, ~sum(is.na(.)))
#merge ethnicity 
dim(DBs)                                        #502411     10
DBs <- merge(DBs,DB_Cov2,by="IID")
names(DBs)
dim(DBs)                                        #502409     12
DBs %>% count(Ethn_b)                           #0 472611, 1 28899, NA 898
#exclude non-Whites
DBs=subset(DBs, Ethn_b==0 & !is.na(Ethn_b)) 
DBs %>% count(Ethn_b)                           #0 472611
DBs=subset(DBs, select= -c(Ethn,Ethn_b))
names(DBs)

#Missing structure: non-Whites
map(DBs, ~sum(!is.na(.)))

#Dementia frequency after excluding non-Whites
print("Frequency of Alzheimer Dementia:")
table(DBs$DemAlz)
print("Frequency of Vascular Dementia:")
table(DBs$DemVas)
print("Frequency of all Dementia:")
table(DBs$Dementia)

#=====Kinship (Relatedness)
kinship = read_table("/data/Wolfson-PNU-dementia/ukb_helper_files/ukb59138_rel_s488264.dat")
summary(kinship$Kinship)
#Min.    1st Qu. Median  Mean    3rd Qu. Max. 
#0.0442  0.0566  0.0699  0.1185  0.2254  0.4999         n=107116
kinship$Kinship_0.125 <- ifelse(kinship$Kinship >= 0.125, c("higher"), c("lower"))
table(kinship$Kinship_0.125)                            #n=73090/34026
#1st degree relatives (kinship coefficient>=0.125)
Exclu_kinship = kinship %>% filter(Kinship>0.125) %>% select(ID1) %>% rename("IID"="ID1") 
#Exclusion/Filtering
DBs <- DBs  %>% mutate(exclu_Kinship = ifelse(IID %in% Exclu_kinship$IID, "1", "0"))
table(DBs$exclu_Kinship)                              #n=442563/30048
#DBs<-filter(DBs, exclu_Kinship==0)
DBs=subset(DBs, exclu_Kinship==0)
table(DBs$exclu_Kinship)                              #n=442563
DBs=subset(DBs, select= -c(exclu_Kinship))
names(DBs)

#Missing structure: non-1st relatives
map(DBs, ~sum(!is.na(.)))

#Dementia frequency after further exclusion of first degree relatives
print("Frequency of Alzheimer Dementia:")
table(DBs$DemAlz)
print("Frequency of Vascular Dementia:")
table(DBs$DemVas)
print("Frequency of all Dementia:")
table(DBs$Dementia)


#=====Poor heterozygosity/missingness
#Read ukb_sqc_v2.txt file
DB_hmid <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_SQC_outliers.rds")

#Merge het.missing.outliers to the outcome dataset
dim(DBs)                                        #442563     5
DBs <- left_join(DBs,DB_hmid,by="IID")
names(DBs)
dim(DBs)                                        #442563     6
DBs %>% count(het.missing.outliers)                           #0 429097, 1 907, NA 12559
#exclude outliers
DBs=subset(DBs, het.missing.outliers==0 | is.na(het.missing.outliers)) 
DBs %>% count(het.missing.outliers)                           #0 429097, NA 12559
DBs=subset(DBs, select= -c(het.missing.outliers, V25))
names(DBs)

#Missing structure: No poor heterozygosity/missingness
map(DBs, ~sum(!is.na(.)))
#Dementia frequency after further exclusion of subjects with poor heterozygosity/missingness
table(DBs$Dementia)


#=====Covariates: PC 
##Read Cov_LF file
DBc <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")
DBcs=subset(DBc, select= c(IID,age,sex,centre,pc1,pc5,FEV1_b0,FVC_b0))
map(DBcs, ~sum(!is.na(.)))

DBcs <- DBcs %>% mutate(PCmis = ifelse(!is.na(pc1) & !is.na(pc5) & !is.na(age),0,1))
table(DBcs$PCmis)
DBcs <- DBcs %>% mutate(LFmis = ifelse(!is.na(FVC_b0),0,1))
table(DBcs$LFmis)
DBcs=subset(DBcs, select= c(IID,PCmis,LFmis))
str(DBcs)

#Merge Covariates&LF missing to the outcome dataset
dim(DBs)                                        #441656     10
#.DBs <- merge(DBs,DB_hmid,by="IID")
DBs <- left_join(DBs,DBcs,by="IID")
names(DBs)
dim(DBs)                                        #441656     12
#exclude those without genetic data
DBs %>% count(PCmis)                            #0 429097, 1 12559
DBs=subset(DBs, PCmis==0) 
DBs %>% count(PCmis)                            #0 429097
DBs=subset(DBs, select= -c(PCmis))
names(DBs)

#Missing structure: No PC missing
map(DBs, ~sum(!is.na(.)))
#Dementia frequency after further exclusion of subjects without PC (genetic data)
DBs %>% count(Dementia)                         #No missing, so no need for a separate file


#=====Lost to follow-up at baseline (R_LTFU_f191)
table(DBs$LTFU_0)                            #Yes 1057, No 428040
DBs <- DBs %>% mutate(LTF_0i=ifelse(LTFU_0=="Yes" & Dementia != 1, 1,0))
table(DBs$LTF_0i)                            #1 1051, 0 428046
DBs=subset(DBs, LTFU_0=="No" | Dementia == 1) 
DBs=subset(DBs, select= -c(LTF_0i))

#Missing structure: Lost to follow-up @Baseline
map(DBs, ~sum(!is.na(.)))
#Dementia frequency after further exclusion of subjects without PC (genetic data)
DBs %>% count(Dementia)    


#=====Subtype dementia
#Alzheimer's 
DBs %>% count(DemAlz)                             #1 2859, 0 425187
DBs$DemAlz[DBs$DemAlz==0 & DBs$Dementia==1] <- NA
DBs %>% count(DemAlz)                             #1 2859, 0 421241, NA 3946

#Vascular Dementia 
DBs %>% count(DemVas)                             #1 1544, 0 426502
DBs$DemVas[DBs$DemVas==0 & DBs$Dementia==1] <- NA
DBs %>% count(DemVas)                             #1 1544, 0 421241, NA 5261


#=====FID must be first
DBs <- DBs %>% relocate(FID, .before = IID)
str(DBs)


#=====Extracting new Tenure variables, YearsWithDiagnosis and g-factors from updated wrangled file
DBu <-readRDS("/data/Wolfson-PNU-dementia/lungdev/gfactor_and_updateddiag_2023/UKBW_59138r670589_Mohammad_New.rds")

DBus = DBu%>% select(contains("EID"),
                    contains("R_CogGSet2d_t1t2t6t7t9t4t5t8t10_Min3TestsnotNA_CognitiveGFactor_2"),
                    contains("R_CogGSet6d_t1t5t10_Min3TestsnotNA_CognitiveGFactor_2"),
                    contains("R_CogGSet7d_t2t6t7t9_Min3TestsnotNA_CognitiveGFactor_2"), 
                    contains("R_ML_C42C240Xf41270f20002_DementiaTenureAdjusted"),
                    contains("ML_C42C240Xf41270_DementiaDiagyrBin0"),
                    contains("R_ML_C42C240Xf41270f20002_Dementia_Incident"),
                    contains("R_ML_C42C240Xf41270f20002_Dementia_Prevalent")) 
DBus <- rename(DBus, c(IID=EID))
DBus <- rename(DBus, c(DementiaTenureA=R_ML_C42C240Xf41270f20002_DementiaTenureAdjusted))
DBus <- rename(DBus, c(Dementia_Incident=R_ML_C42C240Xf41270f20002_Dementia_Incident))
DBus <- rename(DBus, c(Dementia_Prevalent=R_ML_C42C240Xf41270f20002_Dementia_Prevalent))
str(DBus)
saveRDS(DBus, file="/data/Wolfson-PNU-dementia/lungdev/gfactor_and_updateddiag_2023/UKBW_59138r670589_selected.rds")
#selecting tenure
DBuss = DBus%>% select(IID, DementiaTenureA, Dementia_Incident, Dementia_Prevalent)
export(DBuss,"/data/Wolfson-PNU-dementia/lungdev/gfactor_and_updateddiag_2023/DementiaTenureA.dta")
#Merging
DBs <- left_join(DBs,DBuss,by="IID")
names(DBs)


#=====exclude those with dementia at baseline    
DBs$pDementia_0 <- DBs$Dementia_Prevalent
#Incident all-cause dementia
DBs %>% count(pDementia_0)                           #0: 6602 1:205 NA: 421,239
DBs %>% count(Dementia)                              #6805  421241
DBs$iDementia <- DBs$Dementia
DBs %>% count(iDementia)                             #6805  421241
DBs$iDementia[DBs$pDementia_0==1] <- NA
DBs %>% count(iDementia)                             #6601  421240    NA 205
names(DBs)
#Incident Alzheimer's Disease
DBs %>% count(DemAlz)                              #2859  421241    NA 3946
DBs$iDemAlz <- DBs$DemAlz
DBs %>% count(iDemAlz)                             #2859  421241    NA 3946
DBs$iDemAlz[DBs$pDementia_0==1] <- NA
DBs %>% count(iDemAlz)                             #2785  421240    NA 4021
names(DBs)
#Incident vascular dementia
DBs %>% count(DemVas)                              #1544  421241    NA 5261
DBs$iDemVas <- DBs$DemVas
DBs %>% count(iDemVas)                             #1544  421241    NA 5261
DBs$iDemVas[DBs$pDementia_0==1] <- NA
DBs %>% count(iDemVas)                             #1517  421240    NA 5289
names(DBs)


#=====Saving
print("All dementia variables were saved as: ukb_pheno_Dem.rds")
saveRDS(DBs, file="/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_Dem.rds")
print("Selected outcomes were saved as: ukb_pheno_Dem.txt")
write.table(DBs, file = "/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_Dem.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#=====Covariates: LF
#exclude those without LF data              *--> NOTE: Decided NOT to exclude 
DBs %>% count(LFmis)                           #0 322117, 1 105929
#DBs=subset(DBs, LFmis==0) 
DBs %>% count(LFmis)                           #0 322117
#DBs=subset(DBs, select= -c(LFmis))
names(DBs)
map(DBs, ~sum(!is.na(.)))
#Dementia frequency after further exclusion of subjects without LF data (best measure)
table(DBs$Dementia)

#Missing structure: No LF
map(DBs, ~sum(!is.na(.)))
#Dementia frequency after further exclusion of subjects without PC (genetic data)
DBs %>% count(Dementia)   


###============================================Outcome files
##Read extracted file, managed
DBs <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_Dem.rds")

#1.All Cause Dementia file
Dementia = DBs%>% select(contains("FID"), contains("IID"), Dementia)
map(Dementia, ~sum(!is.na(.)))
Dementia <- Dementia[complete.cases(Dementia),]
str(Dementia)
write.table(Dementia, file = "/data/Wolfson-PNU-dementia/lungdev/Dementia.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#2.AD file
DemAlz = DBs%>% select(contains("FID"), contains("IID"), DemAlz)
map(DemAlz, ~sum(!is.na(.)))
DemAlz <- DemAlz[complete.cases(DemAlz),]
str(DemAlz)
write.table(DemAlz, file = "/data/Wolfson-PNU-dementia/lungdev/DemAlz.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#3.Vascular dementia file
DemVas = DBs%>% select(contains("FID"), contains("IID"), DemVas)
map(DemVas, ~sum(!is.na(.)))
DemVas <- DemVas[complete.cases(DemVas),]
str(DemVas)
write.table(DemVas, file = "/data/Wolfson-PNU-dementia/lungdev/DemVas.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

