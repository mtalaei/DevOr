#-----------------------------------------------------------------------------
#Covariates & Lung function file 
#-----------------------------------------------------------------------------

library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")

#managed file
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")


#-----read Covariates
print("-----------------------------------------*** ukb_cov ***")
DB_Cov <-readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_cov.rds")
#DB_Cov <- DB_Cov[-c(1)]         #Omit FID
print("*** Structure- Covariates:")
dim(DB_Cov)
str(DB_Cov)
print("-----------------------------------------ukb_cov ***End")

#-----read Covariates 2
print("-----------------------------------------*** CoreDataClin1 ***")
DB <- read_tsv("/data/Wolfson-PNU-dementia/datasets/ukb_extracted/47532_Coded/1_ukbpheno_CoreDataClin1.tsv")
DB_Cov2 = DB%>% select(contains("EID"),
                   contains("Standing height.0.0"),
                   contains("Ethnic background.0.0"),
                   contains("Townsend deprivation index at recruitment.0.0"),
                   contains("Qualifications.0.0"))
DB_Cov2 <- rename(DB_Cov2, "IID"="EID")
DB_Cov2 <- rename(DB_Cov2, "Height_0"="Standing height.0.0")
DB_Cov2 <- rename(DB_Cov2, "Ethn"="Ethnic background.0.0")
DB_Cov2 <- rename(DB_Cov2, "TDI"="Townsend deprivation index at recruitment.0.0")
DB_Cov2 <- rename(DB_Cov2, "Edu_0"="Qualifications.0.0")
DB_Cov2 <- DB_Cov2[-c(3)]
print("----Selection structure")
str(DB_Cov2)
print("----Selection names")
names(DB_Cov2)
print("----Standing height.0.0")
describe(DB_Cov2$Height_0)
print("----Ethnic background")
table(DB_Cov2$Ethn)
DB_Cov2 %>% count(Ethn)
print("----Townsend deprivation index")
describe(DB_Cov2$TDI)
print("----Qualifications")
table(DB_Cov2$Edu_0)
DB_Cov2 %>% count(Edu_0)
print("-----------------------------------------CoreDataClin1 ***End")


#-----read Smoking
print("-----------------------------------------*** 11_ukbpheno_AlcSmoke ***")
DB <- read_tsv("/data/Wolfson-PNU-dementia/datasets/ukb_extracted/47532_Coded/11_ukbpheno_AlcSmoke.tsv")
DB_Sm = DB%>% select(contains("EID"), 
                    contains("Smoking status.0.0"), 
                    contains("Pack years of smoking.0.0"))
DB_Sm <- rename(DB_Sm, "IID"="EID")
DB_Sm <- rename(DB_Sm, "Smoking_0"="Smoking status.0.0")
DB_Sm <- rename(DB_Sm, "Packyears_0"="Pack years of smoking.0.0")
print("*** Structure- Smoking:")
dim(DB_Sm)
str(DB_Sm)
print("----Smoking status")
table(DB_Sm$Smoking_0)
print("----Pack years of smoking")
describe(DB_Sm$Packyears_0)
print("-----------------------------------------11_ukbpheno_AlcSmoke ***End")


#-----read LF (Best measure at instance 0)
print("-----------------------------------------*** 9_ukbpheno_PhysMsr ***")
DB <- read_tsv("/data/Wolfson-PNU-dementia/datasets/ukb_extracted/47532_Coded/9_ukbpheno_PhysMsr.tsv")
DB_LF = DB%>% select(contains("EID"), contains("Best measure.0.0"))
DB_LF <- rename(DB_LF, "IID"="EID")
DB_LF <- rename(DB_LF, "FEV1_b0"="Forced expiratory volume in 1-second (FEV1), Best measure.0.0")
DB_LF <- rename(DB_LF, "FVC_b0"="Forced vital capacity (FVC), Best measure.0.0")
DB_LF$FEV1FVC_b0 <- DB_LF$FEV1_b0/DB_LF$FVC_b0 
print("*** Structure- Lung function instance 0 (Best):")
dim(DB_LF)
str(DB_LF)
print("----FEV1_best summary")
describe(DB_LF$FEV1_b0)
print("----FVC_best summary")
describe(DB_LF$FVC_b0)
print("----FEV1/FVC_best summary")
describe(DB_LF$FEV1FVC_b0)
print("-----------------------------------------9_ukbpheno_PhysMsr ***End")

#-----Merging
DB <- merge(DB_LF,DB_Cov,by="IID")
DB <- merge(DB,DB_Cov2,by="IID")
DB <- merge(DB,DB_Sm,by="IID")
print("*** Structure- Merged dataset:")
dim(DB)
str(DB)
print("*** Selection names")
names(DB)


#-----adding age^2 and moving FID
DB <- DB %>% mutate(age2=age*age)
DB <- DB %>% relocate(FID, .before = IID)
str(DB)





#-----Saving
saveRDS(DB, file="/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")
#selection (for covariate file to be used by Regenie)
DBs = DB%>% subset(select=c(FID, IID, age, age2, sex, array, centre, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10, 
                            Height_0, Smoking_0))
map(DBs, ~sum(!is.na(.)))
DBs <- DBs[complete.cases(DBs),]
write.table(DBs, file = "/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#-----------------------------------------------------------------------------
#Lancet Commission variables
#-----------------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) 


library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")

#=====read Sheena's extracted file
##Read UKBW_78867r675516_v290524_Mohammad.rds
DB_C <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/UKBW_78867r675516_v100624_Mohammad.rds")
#Explore
DBs = DB_C%>% select(contains("EID"),
                   contains("20256"),
                   contains("20258"),
                   contains("f20022"),
                   contains("DementiaTenure"))
names(DBs)
str(DBs)
describe(DBs$R_ML_C42C240Xf41270f20002_DementiaTenureAdjusted)
                
##select of 34 covarates (best version)
DB_Cs=subset(DB_C, select= c(EID59138, EID78867,
                             R_RF_DemLC_1_LackSecondaryEdu,
                             R_RF_DemLC_2_Hearing_0,
                             R_RF_DemLC_3_HeadInjury_0,
                             R_RF_DemLC_4_HypertensionMidLife_0,
                             R_RF_DemLC_5_Alcohol_0,
                             R_RF_DemLC_6_BMI_obese_0,
                             R_RF_DemLC_7_Diabetes_0,
                             R_RF_DemLC_8_SmokingCurrent_0,
                             R_RF_DemLC_9_SocialIsolation_0,
                             R_RF_DemLC_10_Depression_0,
                             R_RF_DemLC_11_PhysicalInactivity_0,
                             R_RF_DemLC_12_AirPollution,
                             APOE4_alleles,
                             R_BMI_f21001_0,
                             ipaq_activity_group_f22032_0_0,
                             ML_C42C240Xf41270f20002_IschaemicHeart, 
                             R_ML_DiagbySess_C42C240Xf41270f20002_IschaemicHeart_0,
                             ML_C42C240Xf41270f20002_HeartFailure, 
                             R_ML_DiagbySess_C42C240Xf41270f20002_HeartFailure_0,
                             ML_C42C240Xf41270f20002_Stroke,
                             R_ML_DiagbySess_C42C240Xf41270f20002_Stroke_0,
                             ML_C42C240Xf41270f20002_Asthma,
                             R_ML_DiagbySess_C42C240Xf41270f20002_Asthma_0,
                             ML_C42C240Xf41270_BreathingProblems,
                             ML_C42C240Xf41270f20002_EmphysemaChronBronchitis,
                             ML_C2410_J47bronchiectasis,
                             R_ML_CholesterolZComposite_f6153f6177f30690C42C240Xf41270f20002_0,
                             forced_vital_capacity_fvc_zscore_f20257_0_0,
                             forced_expiratory_volume_in_1second_fev1_zscore_f20256_0_0,
                             fev1_fvc_ratio_zscore_f20258_0_0,
                             birth_weight_f20022_0_0,
                             birth_weight_f20022_2_0,
                             R_BirthWeight_f20022_0,
                             R_BirthWeight_f20022_2))

dim(DB_Cs)
str(DB_Cs)

#-----Saving
#R
saveRDS(DB_Cs, file="/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_Cov02_LF.rds")

#Stata
DB_Cs <- rename(DB_Cs, c(LC1_no2ndEdu="R_RF_DemLC_1_LackSecondaryEdu"))
DB_Cs <- rename(DB_Cs, c(LC2_Hearing_0="R_RF_DemLC_2_Hearing_0"))
DB_Cs <- rename(DB_Cs, c(LC3_HeadInjury_0="R_RF_DemLC_3_HeadInjury_0"))
DB_Cs <- rename(DB_Cs, c(LC4_HTNmLife_0="R_RF_DemLC_4_HypertensionMidLife_0"))
DB_Cs <- rename(DB_Cs, c(LC5_Alcohol_0="R_RF_DemLC_5_Alcohol_0"))
DB_Cs <- rename(DB_Cs, c(LC6_Obese_0="R_RF_DemLC_6_BMI_obese_0"))
DB_Cs <- rename(DB_Cs, c(LC7_T2D_0="R_RF_DemLC_7_Diabetes_0"))
DB_Cs <- rename(DB_Cs, c(LC8_SmokingC_0="R_RF_DemLC_8_SmokingCurrent_0"))
DB_Cs <- rename(DB_Cs, c(LC9_SocialIso_0="R_RF_DemLC_9_SocialIsolation_0"))
DB_Cs <- rename(DB_Cs, c(LC10_Depression_0="R_RF_DemLC_10_Depression_0"))
DB_Cs <- rename(DB_Cs, c(LC11_PhInactivity_0="R_RF_DemLC_11_PhysicalInactivity_0"))
DB_Cs <- rename(DB_Cs, c(LC12_AirPollution="R_RF_DemLC_12_AirPollution"))
DB_Cs <- rename(DB_Cs, c(BMI_0="R_BMI_f21001_0"))
DB_Cs <- rename(DB_Cs, c(Activity_0="ipaq_activity_group_f22032_0_0"))
DB_Cs <- rename(DB_Cs, c(IschaemicHeart="ML_C42C240Xf41270f20002_IschaemicHeart"))
DB_Cs <- rename(DB_Cs, c(IschaemicHeart_0="R_ML_DiagbySess_C42C240Xf41270f20002_IschaemicHeart_0"))
DB_Cs <- rename(DB_Cs, c(HeartFailure="ML_C42C240Xf41270f20002_HeartFailure"))
DB_Cs <- rename(DB_Cs, c(HeartFailure_0="R_ML_DiagbySess_C42C240Xf41270f20002_HeartFailure_0"))
DB_Cs <- rename(DB_Cs, c(Stroke="ML_C42C240Xf41270f20002_Stroke"))
DB_Cs <- rename(DB_Cs, c(Stroke_0="R_ML_DiagbySess_C42C240Xf41270f20002_Stroke_0"))
DB_Cs <- rename(DB_Cs, c(Asthma="ML_C42C240Xf41270f20002_Asthma"))
DB_Cs <- rename(DB_Cs, c(Asthma_0="R_ML_DiagbySess_C42C240Xf41270f20002_Asthma_0"))
DB_Cs <- rename(DB_Cs, c(BreathingProblems="ML_C42C240Xf41270_BreathingProblems"))
DB_Cs <- rename(DB_Cs, c(EmphysemaChronBronchitis="ML_C42C240Xf41270f20002_EmphysemaChronBronchitis"))
DB_Cs <- rename(DB_Cs, c(Bronchiectasis="ML_C2410_J47bronchiectasis"))
DB_Cs <- rename(DB_Cs, c(CholesterolZComposite_0="R_ML_CholesterolZComposite_f6153f6177f30690C42C240Xf41270f20002_0"))
DB_Cs <- rename(DB_Cs, c(FVCz_0="forced_vital_capacity_fvc_zscore_f20257_0_0"))
DB_Cs <- rename(DB_Cs, c(FEVz_0="forced_expiratory_volume_in_1second_fev1_zscore_f20256_0_0"))
DB_Cs <- rename(DB_Cs, c(FERz_0="fev1_fvc_ratio_zscore_f20258_0_0"))
DB_Cs <- rename(DB_Cs, c(Bweight_0="R_BirthWeight_f20022_0"))
DB_Cs <- rename(DB_Cs, c(Bweight_2="R_BirthWeight_f20022_2"))
DB_Cs <- rename(DB_Cs, c(rBweight_0="birth_weight_f20022_0_0"))
DB_Cs <- rename(DB_Cs, c(rBweight_2="birth_weight_f20022_2_0"))

str(DB_Cs)

export(DB_Cs,"/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_Cov02_LF.dta")


