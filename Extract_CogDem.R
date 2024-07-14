#------------------------------------------------------------------------------
##Extracting Dementia and cognition tests plus HTN, from Wrangled data by Sheena 
#------------------------------------------------------------------------------
library(readr)
library(dplyr)

DB <- readRDS("/data/Wolfson-PNU-dementia/UKB/datasets/sheenastrux/59138r52532/Wrangling/52532_080423/UKBW_59138r52532_S21A.rds")
print("Dimension of Wrangled dataset")
dim(DB)
DBs = DB%>% select(contains("EID"),
                   contains("Dementia"),
                   contains("Alzheimers"),
                   contains("20023"), 
                   contains("4282"),
                   contains("20016"), 
                   contains("399"),
                   contains("20018"),
                   contains("6348"),
                   contains("6350"),
                   contains("23324"),
                   contains("6373"),
                   contains("21004"),
                   contains("20197"),
                   contains("ML_C42C240Xf41270_DementiaDiagyrBin"),
                   contains("ML_C42C240Xf41270_DementiaAge"),
                   contains("ML_C42C240Xf41270_DementiaDiagyr"),
                   contains("R_Age_f21003"),
                   contains("R_Date_f53"),
                   contains("R_LTFU"),
                   contains("R_Death"),
                   contains("R_ExclusionWithdrawals"),
                   contains("4079"),
                   contains("4080"),
                   contains("6177"),
                   contains("ML_f41270_I10Hypertension"),
                   contains("ML_f20002_SR1072EssentialHypertension"),
                   contains("C2409_I10essentialprimaryHypertension"),
                   contains("ML_C2409_I10essentialprimaryHypertensionAge"),
                   contains("ML_f41270_I10HypertensionAge"),
                   contains("ML_f20002_SR1072EssentialHypertensionAge"),
                   contains("ML_C42C240Xf41270f20002_Asthma"), 
                   contains("ML_C42C240Xf41270f20002_IschaemicHeartMI"),
                   contains("ML_C42C240Xf41270f20002_Angina"),
                   contains("ML_f20004_OP1069Heart"),
                   contains("ML_f20004_OP1070Angioplasty"),
                   contains("ML_f20004_OP1095CoronaryArteryBypass"),
                   contains("R_ML_IschaemicHeartMIComposite_f6150f3894C42C240Xf41270f20002"))
print("Dimension of cognition/dementia selected dataset")
str(DBs)
dim(DBs)
#Saving
print("All cognition, dementia, and hypertension variables were saved as: UKBW_52532_110822A_CogDem.rds")
saveRDS(DBs, file="/data/Wolfson-PNU-dementia/lungdev/UKBW_52532_110822A_CogDem.rds")
