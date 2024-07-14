#-----------------------------------------------------------------------------
#Lung function - Cognitive function/Dementia file for testing associations
#-----------------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) 


library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")

##Set working directory 
setwd("/data/home/hmy431/scripts")

#=====read CF (selected cognitive function variables at instance 2) --ukb_pheno_cogT2S
DB_CF <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_cogTS.rds")

DB_CF=subset(DB_CF, select= c(IID,
                              PairsMatching_0_logrz,	PairsMatching_2_logrz, 
                              ReactionTime_0_logrz,	ReactionTime_2_logrz,
                              ProspectiveMemory_0_b,	ProspectiveMemory_2_b, 
                              FluidIntelligence_0_z,	FluidIntelligence_2_z,
                              NumericMemory_0_z,	NumericMemory_2_z, 
                              TrailMakingN_du_2_logrz,
                              TrailMakingAN_du_2_logrz, 
                              TrailMakingANN_du_2_rz,
                              TrailMakingANNr_du_2_logrz, 
                              SymbolDigit_c_2_z,
                              MatrixPattern_2_z,
                              TowerRearranging_c_2_z, 
                              PairedLearning_2_z,
                              gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
print("*** Structure- Cognitive function instance 2:")
dim(DB_CF)
str(DB_CF)
print("-----------------------------------------ukb_pheno_cogTS ***End")

#=====read Dm (Dementia and its subtypes)
DB_Dm <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_Dem.rds")
DB_Dm <- subset(DB_Dm,select=-c(LFmis, FID))           #Omit LFmis
str(DB_Dm)
DB_Dm %>% count(Dementia) 


#=====read Cov_LF: All covariates (age,sex,centre,ethnicity,height,Townsend deprivation index,smoking,PC & LF
##Read Cov_LF file
DB_C <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")
DB_C <- subset(DB_C,select=-c(FID))           #Omit FID
str(DB_C)

map(DB_C, ~sum(!is.na(.)))
#DB_Cov <- DB_Cov[-c(1)]         #Omit FID


##=====HTN 
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/UKBW_52532_110822A_CogDem.rds")
##Select components
DB_H = DB%>% select("EID", 
                   "R_Systolic_f4080_mean_0",
                   "systolic_blood_pressure_automated_reading_f4080_0_1",
                   "R_Diastolic_f4079_mean_0",
                   "diastolic_blood_pressure_automated_reading_f4079_0_1",
                   "R_MedBloodPressureComposite_f6153f6177_0",
                   "ML_f41270_I10Hypertension", 
                   "ML_f20002_SR1072EssentialHypertension", 
                   "ML_C2409_I10essentialprimaryHypertension",
                   "ML_C2409_I10essentialprimaryHypertensionAge",
                   "R_Age_f21003_0")
names(DB_H)

#Renaming
DB_H <- rename(DB_H, c(IID=EID))
DB_H$FID <- DB_H$IID
DB_H <- DB_H %>% relocate(FID, .before = IID)
str(DB_H)

#Defining HTN at baseline
#Note: mutate(R_Systolic_f4080_mean_0=rowMeans(cbind(systolic_blood_pressure_automated_reading_f4080_0_0,systolic_blood_pressure_automated_reading_f4080_0_1), na.rm = TRUE))
describe(DB_H$R_Systolic_f4080_mean_0)                          #min 65 max 268  n 472294
describe(DB_H$R_Diastolic_f4079_mean_0)                         #min 32 max 147.5  n 472299
table(DB_H$R_MedBloodPressureComposite_f6153f6177_0)            #No 389,819    Yes 103,987 Unknown 8605
table(DB_H$ML_C2409_I10essentialprimaryHypertension)            #0 307460    1 194951
#mean BP and 6177 + 6153
DB_H <- DB_H%>% mutate(HTN_02=ifelse(R_Diastolic_f4079_mean_0 >= 90 |
                                       R_Systolic_f4080_mean_0 >= 140 |
                                       R_MedBloodPressureComposite_f6153f6177_0 =="Yes",1,0))
table(DB_H$HTN_02)                                #0 220,454 1 257,781 
#2nd BP and 6177 + 6153
DB_H <- DB_H%>% mutate(HTN_03=ifelse(diastolic_blood_pressure_automated_reading_f4079_0_1 >= 90 |
                                       systolic_blood_pressure_automated_reading_f4080_0_1 >= 140 |
                                       R_MedBloodPressureComposite_f6153f6177_0 =="Yes",1,0))
table(DB_H$HTN_03)                                #0 222,184  1 247,438 
#2nd BP and ICD10 (damaged data)
DB_H <- DB_H %>% mutate(pHTN_I10_0 = ifelse(ML_C2409_I10essentialprimaryHypertensionAge<=R_Age_f21003_0,1,0))
table(DB_H$pHTN_I10_0)
DB_H <- DB_H%>% mutate(HTN_04=ifelse(diastolic_blood_pressure_automated_reading_f4079_0_1 >= 90 |
                                       systolic_blood_pressure_automated_reading_f4080_0_1 >= 140 |
                                       pHTN_I10_0 ==1,1,0))
table(DB_H$HTN_04) 
table(DB_H$HTN_02, DB_H$HTN_03)

#missing structure
map(DB_H, ~sum(!is.na(.)))

#FO for Essential Hypertension            ---> NOT suitable to use at baseline (covers the whole period)
table(DB_H$HTN_0, DB_H$ML_C2409_I10essentialprimaryHypertension) 

DB_Hs = DB_H%>% select("IID", "HTN_03")
DB_Hs = DB_Hs%>% rename(HTN_0="HTN_03")
str(DB_Hs)
table(DB_Hs$HTN_0)                                #0 222,184  1 247,438


##=====Merging
DBm <- full_join(DB_Dm,DB_CF,by="IID")
DBm <- full_join(DBm,DB_C,by="IID")
DBm <- full_join(DBm,DB_Hs,by="IID")

print("*** Structure- Merged dataset:")
dim(DBm)                                          #502464     57
str(DBm)
print("*** Selection names")
names(DBm)


##=====Incident Alzheimer's Disease and Vascular dementia
#Note: The output file does not include iDemAlz.txt iDemVas.txt because they were created later (No need in case of repeating from scratch) 
DBm %>% count(DemAlz)                              #2859  421241    NA 78364
DBm$iDemAlz <- DBm$DemAlz
DBm %>% count(iDemAlz)                             #2859  421241    NA 78364
DBm$iDemAlz[DBm$pDementia_0==1] <- NA
DBm %>% count(iDemAlz)                             #2841  421241    NA 78382
names(DBm)
#Incident vascular dementia
DBm %>% count(DemVas)                              #1544  421241    NA 79679
DBm$iDemVas <- DBm$DemVas
DBm %>% count(iDemVas)                             #1544  421241    NA 79679
DBm$iDemVas[DBm$pDementia_0==1] <- NA
DBm %>% count(iDemVas)                             #1533  421241    NA 79690
names(DBm)


##Modified Tenure variable was added to ukb_pheno_Dem.rds and ukb_LF0-DmCF_m.dta


#-----Saving
saveRDS(DBm, file="/data/Wolfson-PNU-dementia/lungdev/ukb_LF0-DmCF.rds")
export(DBm,"/data/Wolfson-PNU-dementia/lungdev/ukb_LF0-DmCF.dta")





#-----------------------------------------------------------------------------
#Lancet Commission variables
#-----------------------------------------------------------------------------
#=====read Sheena's extracted file
##Read UKBW_78867r675516_v290524_Mohammad.rds
DB_C <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_Cov02_LF.rds")

##Note: variables extracted in Cov_LF_file.R and was merged in Stata: DaMa_LF-CF.do



