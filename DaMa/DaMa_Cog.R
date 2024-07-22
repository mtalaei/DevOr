#-------------------------------------------------------------------------------
# Cognition data management (instance 0 @Baseline AND 2)
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) 

library(tidyverse)
library(readr)
library(dplyr)
library(psych)
library("rio")
library(naniar)
library(ukbtools)


##====Read 
#raw extracted file
DB <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/UKBW_52532_110822A_CogDem.rds")

#managed file
DBs <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_pheno_cogTS.rds")

#Explore: 399, 20023, 20018, 20016, 4282--------------------------------------- 
DBs = DB%>% select(contains("4282"))
names(DBs)
hist(DBs$R_Test3_ProspectiveMemory_f20018_0, breaks=30)
DBs%>%count(ML_f41270_I10Hypertension)
#------------------------------------------------------------------------------

##====Selected variables at instances 0 and 2
DBs = DB%>% select("EID", 
                   contains("R_Test1_PairsMatching_f399_0"),
                   contains("R_Test1_PairsMatching_f399_2"),
                   contains("R_Test2_ReactionTime_f20023_0"),
                   contains("R_Test2_ReactionTime_f20023_2"), 
                   contains("R_Test3_ProspectiveMemory_f20018_0"), 
                   contains("R_Test3_ProspectiveMemory_f20018_2"), 
                   contains("R_Test4_FluidIntelligence_f20016_0"),  
                   contains("R_Test4_FluidIntelligence_f20016_2"),  
                   contains("R_Test5_NumericMemory_f4282_0"),
                   contains("R_Test5_NumericMemory_f4282_2"),
                   contains("R_Test6_TrailMaking_ANRT_f6350_2"),
                   contains("R_Test6_TrailMaking_NRT_f6348_2"),
                   contains("R_Test7_DigitSymbol_NCFinal_f23324_2"),
                   contains("R_Test8_MatrixPattern_NCFinal_f6373_2"),
                   contains("R_Test8_MatrixPattern_NCFinalNARat_f6373f6374_2"),
                   contains("R_Test9_TowerRearranging_NCFinal_f21004_2"),
                   contains("R_Test10_PairedAssociate_f20197_2"),
                   contains("ML_C42C240Xf41270_DementiaDiagyrBin0"),
                   contains("ML_C42C240Xf41270_DementiaAge"),
                   contains("R_Age_f21003_2"),
                   contains("R_ExclusionWithdrawals"),
                   contains("R_LTFU_f191"))
names(DBs)

#Renaming
DBs <- rename(DBs, c(IID=EID))
DBs$FID <- DBs$IID
DBs <- DBs %>% relocate(FID, .before = IID)
print("Structure of selected CF variables @Baseline:")
str(DBs)

DBs <- rename(DBs, 
              PairsMatching_0='R_Test1_PairsMatching_f399_0',
              PairsMatching_2='R_Test1_PairsMatching_f399_2',              
              ReactionTime_0= 'R_Test2_ReactionTime_f20023_0',
              ReactionTime_2= 'R_Test2_ReactionTime_f20023_2',              
              ProspectiveMemory_0= 'R_Test3_ProspectiveMemory_f20018_0',
              ProspectiveMemory_2= 'R_Test3_ProspectiveMemory_f20018_2',              
              FluidIntelligence_0= 'R_Test4_FluidIntelligence_f20016_0',
              FluidIntelligence_2= 'R_Test4_FluidIntelligence_f20016_2',              
              NumericMemory_0= 'R_Test5_NumericMemory_f4282_0',
              NumericMemory_2= 'R_Test5_NumericMemory_f4282_2',
              TrailMakingAN_du_2= 'R_Test6_TrailMaking_ANRT_f6350_2',
              TrailMakingN_du_2= 'R_Test6_TrailMaking_NRT_f6348_2',
              SymbolDigit_c_2= 'R_Test7_DigitSymbol_NCFinal_f23324_2',
              MatrixPattern_c_2= 'R_Test8_MatrixPattern_NCFinal_f6373_2',
              MatrixPattern_2= 'R_Test8_MatrixPattern_NCFinalNARat_f6373f6374_2',
              TowerRearranging_c_2= 'R_Test9_TowerRearranging_NCFinal_f21004_2',
              PairedLearning_2= 'R_Test10_PairedAssociate_f20197_2')
str(DBs)              #502411 obs. of  24 variables

##Set working directory for outputs (figures)
setwd("/data/home/hmy431/outputs/Figures/")


#1.1====PairsMatching_0
table(DBs$PairsMatching_0) 
#Log transforming
DBs <- DBs %>% mutate(PairsMatching_0_log=log(PairsMatching_0+1))
describe(DBs$PairsMatching_0_log)           #min 0 max 4.99
#Reverse log
DBs <- DBs %>% mutate(PairsMatching_0_logr=4.99-PairsMatching_0_log)
#Z score reverse log
DBs <- DBs %>% mutate(PairsMatching_0_logrz = (PairsMatching_0_logr - mean(PairsMatching_0_logr, na.rm = TRUE))/sd(PairsMatching_0_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(PairsMatching_0_st=ifelse(!is.na(PairsMatching_0),1,0))

#Report
describe(DBs$PairsMatching_0)               #min=0, max=146
describe(DBs$PairsMatching_0_log)
describe(DBs$PairsMatching_0_logr)
describe(DBs$PairsMatching_0_logrz)
#table(DBs$PairsMatching_0_st) 

#Figures
jpeg(filename = "1_PairsMatching_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_0, breaks=20)          #skewed
dev.off() 
jpeg(filename = "1_PairsMatching_0_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_0_log, breaks=20)      #OK, dentate
dev.off()
jpeg(filename = "1_PairsMatching_0_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_0_logrz, breaks=20)
dev.off() 

#1.2====PairsMatching_2
table(DBs$PairsMatching_2) 
#Log transforming
DBs <- DBs %>% mutate(PairsMatching_2_log=log(PairsMatching_2+1))
describe(DBs$PairsMatching_2_log)           #min 0 max 3.81
#Reverse log
DBs <- DBs %>% mutate(PairsMatching_2_logr=3.81-PairsMatching_2_log)
#Z score reverse log
DBs <- DBs %>% mutate(PairsMatching_2_logrz = (PairsMatching_2_logr - mean(PairsMatching_2_logr, na.rm = TRUE))/sd(PairsMatching_2_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(PairsMatching_2_st=ifelse(!is.na(PairsMatching_2),1,0))

#Report
describe(DBs$PairsMatching_2)               #min=0, max=44
describe(DBs$PairsMatching_2_log)
describe(DBs$PairsMatching_2_logr)
describe(DBs$PairsMatching_2_logrz)
#table(DBs$PairsMatching_2_st) 

#Figures
jpeg(filename = "1_PairsMatching_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_2, breaks=20)          #skewed
dev.off() 
jpeg(filename = "1_PairsMatching_2_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_2_log, breaks=20)      #OK, detate
dev.off()
jpeg(filename = "1_PairsMatching_2_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairsMatching_2_logrz, breaks=20)
dev.off() 


#2.1====ReactionTime_0  
table(DBs$ReactionTime_0)   #there are none under 100 and none over 2000  
#Log transforming
DBs <- DBs %>% mutate(ReactionTime_0_log=log(ReactionTime_0))
describe(DBs$ReactionTime_0_log)              #min=4.62, max=7.6
#Reverse log
DBs <- DBs %>% mutate(ReactionTime_0_logr=7.6-ReactionTime_0_log)
#Z score reverse log
DBs <- DBs %>% mutate(ReactionTime_0_logrz = (ReactionTime_0_logr - mean(ReactionTime_0_logr, na.rm = TRUE))/sd(ReactionTime_0_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(ReactionTime_0_st=ifelse(!is.na(ReactionTime_0),1,0))

#Report
describe(DBs$ReactionTime_0)                #min=101, max=1999
describe(DBs$ReactionTime_0_log)
describe(DBs$ReactionTime_0_logr)
describe(DBs$ReactionTime_0_logrz)
#table(DBs$ReactionTime_0_st) 

#Figures
jpeg(filename = "2_ReactionTime_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_0, breaks=20)         #skewed
dev.off() 
jpeg(filename = "2_ReactionTime_0_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_0_log, breaks=20)
dev.off()
jpeg(filename = "2_ReactionTime_0_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_0_logrz, breaks=20)
dev.off() 

#2.2====ReactionTime_2  
table(DBs$ReactionTime_2)   #there are none under 100 and none over 2000 at instance_2 
#Log transforming
DBs <- DBs %>% mutate(ReactionTime_2_log=log(ReactionTime_2))
describe(DBs$ReactionTime_2_log)              #min=5.01, max=7.5
#Reverse log
DBs <- DBs %>% mutate(ReactionTime_2_logr=7.5-ReactionTime_2_log)
#Z score reverse log
DBs <- DBs %>% mutate(ReactionTime_2_logrz = (ReactionTime_2_logr - mean(ReactionTime_2_logr, na.rm = TRUE))/sd(ReactionTime_2_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(ReactionTime_2_st=ifelse(!is.na(ReactionTime_2),1,0))

#Report
describe(DBs$ReactionTime_2)                #min=150, max=1809
describe(DBs$ReactionTime_2_log)
describe(DBs$ReactionTime_2_logr)
describe(DBs$ReactionTime_2_logrz)
#table(DBs$ReactionTime_2_st) 

#Figures
jpeg(filename = "2_ReactionTime_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_2, breaks=20)         #skewed
dev.off() 
jpeg(filename = "2_ReactionTime_2_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_2_log, breaks=20)
dev.off()
jpeg(filename = "2_ReactionTime_2_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$ReactionTime_2_logrz, breaks=20)
dev.off() 


#3.1====ProspectiveMemory_0
DBs <- DBs %>% mutate(ProspectiveMemory_0_b = recode(ProspectiveMemory_0, 
                                                     "Instruction not recalled, either skipped or incorrect"=1, 
                                                     "Correct recall on second attempt"=1, 
                                                     "Correct recall on first attempt"=0))
table(DBs$ProspectiveMemory_0_b)                #(0)130858  (1)40656
table(DBs$ProspectiveMemory_0, DBs$ProspectiveMemory_0_b)
#Reverse
DBs <- DBs %>% mutate(ProspectiveMemory_0_br = recode(ProspectiveMemory_0_b, "1"="0", "0"="1"))
table(DBs$ProspectiveMemory_0, DBs$ProspectiveMemory_0_br)


#3.2====ProspectiveMemory_2
table(DBs$ProspectiveMemory_2)
DBs <- DBs %>% mutate(ProspectiveMemory_2_b = recode(ProspectiveMemory_2, 
                                                     "Instruction not recalled, either skipped or incorrect"=1, 
                                                     "Correct recall on second attempt"=1, 
                                                     "Correct recall on first attempt"=0))
table(DBs$ProspectiveMemory_2_b)                #(0)41745  (1)8977
table(DBs$ProspectiveMemory_2, DBs$ProspectiveMemory_2_b)
#Reverse
DBs <- DBs %>% mutate(ProspectiveMemory_2_br = recode(ProspectiveMemory_2_b, "1"="0", "0"="1"))
table(DBs$ProspectiveMemory_2, DBs$ProspectiveMemory_2_br)


#4.1====FluidIntelligence_0
table(DBs$FluidIntelligence_0) 
#Z score
DBs <- DBs %>% mutate(FluidIntelligence_0_z = (FluidIntelligence_0 - mean(FluidIntelligence_0, na.rm = TRUE))/sd(FluidIntelligence_0, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(FluidIntelligence_0_st=ifelse(!is.na(FluidIntelligence_0),1,0))

#Report
describe(DBs$FluidIntelligence_0)                #min=0, max=13
describe(DBs$FluidIntelligence_0_z)
#table(DBs$FluidIntelligence_0_st)

#Figures
jpeg(filename = "4_FluidIntelligence_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$FluidIntelligence_0, breaks=20)         #beautifully normal
dev.off() 
jpeg(filename = "4_FluidIntelligence_0_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$FluidIntelligence_0_z, breaks=20)         
dev.off() 


#4.2====FluidIntelligence_2
table(DBs$FluidIntelligence_2) 
#Z score
DBs <- DBs %>% mutate(FluidIntelligence_2_z = (FluidIntelligence_2 - mean(FluidIntelligence_2, na.rm = TRUE))/sd(FluidIntelligence_2, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(FluidIntelligence_2_st=ifelse(!is.na(FluidIntelligence_2),1,0))

#Report
describe(DBs$FluidIntelligence_2)                #min=0, max=13
describe(DBs$FluidIntelligence_2_z)
#table(DBs$FluidIntelligence_2_st)

#Figures
jpeg(filename = "4_FluidIntelligence_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$FluidIntelligence_2, breaks=20)         #beautifully normal
dev.off() 
jpeg(filename = "4_FluidIntelligence_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$FluidIntelligence_2_z, breaks=20)         
dev.off() 


#5.1====NumericMemory_0
table(DBs$NumericMemory_0) 
hist(DBs$NumericMemory_0, breaks=20)
#Z score
DBs <- DBs %>% mutate(NumericMemory_0_z = (NumericMemory_0 - mean(NumericMemory_0, na.rm = TRUE))/sd(NumericMemory_0, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(NumericMemory_0_st=ifelse(!is.na(NumericMemory_0),1,0))

#Report
describe(DBs$NumericMemory_0)                #min=2, max=12
describe(DBs$NumericMemory_0_z)
#table(DBs$NumericMemory_0_st)

#Figures
jpeg(filename = "5_NumericMemory_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$NumericMemory_0, breaks=20)         #beautifully normal
dev.off() 
jpeg(filename = "5_NumericMemory_0_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$NumericMemory_0_z, breaks=20)         
dev.off() 


#5.2====NumericMemory_2
table(DBs$NumericMemory_2) 
hist(DBs$NumericMemory_2, breaks=20)
#Z score
DBs <- DBs %>% mutate(NumericMemory_2_z = (NumericMemory_2 - mean(NumericMemory_2, na.rm = TRUE))/sd(NumericMemory_2, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(NumericMemory_2_st=ifelse(!is.na(NumericMemory_2),1,0))

#Report
describe(DBs$NumericMemory_2)                #min=2, max=12
describe(DBs$NumericMemory_2_z)
#table(DBs$NumericMemory_2_st)

#Figures
jpeg(filename = "5_NumericMemory_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$NumericMemory_2, breaks=20)         #beautifully normal
dev.off() 
jpeg(filename = "5_NumericMemory_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$NumericMemory_2_z, breaks=20)         
dev.off() 


#6====TrailMaking           duration_to_complete_alphanumeric_path_trail_2 ?
table(DBs$TrailMakingN_du_2)
table(DBs$TrailMakingAN_du_2) 

#6.1==TrailMakingN_du_2        A: Numeric
#Log transforming
DBs <- DBs %>% mutate(TrailMakingN_du_2_log=log(TrailMakingN_du_2))
describe(DBs$TrailMakingN_du_2_log)           #min 4.54 max 7.76
#Reverse log
DBs <- DBs %>% mutate(TrailMakingN_du_2_logr=7.76-TrailMakingN_du_2_log)
#Z score reverse log
DBs <- DBs %>% mutate(TrailMakingN_du_2_logrz = (TrailMakingN_du_2_logr - mean(TrailMakingN_du_2_logr, na.rm = TRUE))/sd(TrailMakingN_du_2_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(TrailMakingN_du_2_st = ifelse(!is.na(TrailMakingN_du_2),1,0))

#Report
describe(DBs$TrailMakingN_du_2)                #min=94, max=2356
describe(DBs$TrailMakingN_du_2_log)
describe(DBs$TrailMakingN_du_2_logr)
describe(DBs$TrailMakingN_du_2_logrz)
#table(DBs$TrailMakingN_du_2_st) 

#Figures
jpeg(filename = "61_TrailMakingN_du_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingN_du_2, breaks=20)         #skewed
dev.off() 
jpeg(filename = "61_TrailMakingN_du_2_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingN_du_2_log, breaks=20)
dev.off()
jpeg(filename = "61_TrailMakingN_du_2_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingN_du_2_logrz, breaks=20)
dev.off() 

#6.2==TrailMakingAN_du_2       B: Alphanumeric
#Log transforming
DBs <- DBs %>% mutate(TrailMakingAN_du_2_log=log(TrailMakingAN_du_2))
describe(DBs$TrailMakingAN_du_2_log)           #min 5.08 max 8.66
#Reverse log
DBs <- DBs %>% mutate(TrailMakingAN_du_2_logr=8.66-TrailMakingAN_du_2_log)
#Z score reverse log
DBs <- DBs %>% mutate(TrailMakingAN_du_2_logrz = (TrailMakingAN_du_2_logr - mean(TrailMakingAN_du_2_logr, na.rm = TRUE))/sd(TrailMakingAN_du_2_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(TrailMakingAN_du_2_st = ifelse(!is.na(TrailMakingAN_du_2),1,0))

#Report
describe(DBs$TrailMakingAN_du_2)                #min=160, max=5768
describe(DBs$TrailMakingAN_du_2_log)
describe(DBs$TrailMakingAN_du_2_logr)
describe(DBs$TrailMakingAN_du_2_logrz)
#table(DBs$TrailMakingAN_du_2_st) 

#Figures
jpeg(filename = "62_TrailMakingAN_du_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingAN_du_2, breaks=20)         #skewed
dev.off() 
jpeg(filename = "62_TrailMakingAN_du_2_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingAN_du_2_log, breaks=20)
dev.off()
jpeg(filename = "62_TrailMakingAN_du_2_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingAN_du_2_logrz, breaks=20)
dev.off() 

#6.3==TMT B-A
DBs <- DBs %>% mutate(TrailMakingANN_du_2 = TrailMakingAN_du_2 - TrailMakingN_du_2)
describe(DBs$TrailMakingANN_du_2)                #min=-763, max=5313

#Note: Log transformation is not justified as not much skewness is observed, but there is much kurtosis
#Reverse
DBs <- DBs %>% mutate(TrailMakingANN_du_2_r=5313-TrailMakingANN_du_2)
#Z score reverse log
DBs <- DBs %>% mutate(TrailMakingANN_du_2_rz = (TrailMakingANN_du_2_r - mean(TrailMakingANN_du_2_r, na.rm = TRUE))/sd(TrailMakingANN_du_2_r, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(TrailMakingANN_du_2_st = ifelse(!is.na(TrailMakingANN_du_2),1,0))

#Report
describe(DBs$TrailMakingANN_du_2)                #min -763 max 5313
describe(DBs$TrailMakingANN_du_2_r)
describe(DBs$TrailMakingANN_du_2_rz)
#table(DBs$TrailMakingANN_du_2_st) 

#Figures
jpeg(filename = "63_TrailMakingANN_du_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingANN_du_2, breaks=50)         #skewed
dev.off() 
jpeg(filename = "63_TrailMakingANN_du_2_rz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingANN_du_2_rz, breaks=50)
dev.off() 

#6.4==TMT B/A (AN:N or B:A durations only--measure of executive function)
DBs <- DBs %>% mutate(TrailMakingANNr_du_2 = TrailMakingAN_du_2 / TrailMakingN_du_2)

#Log transforming
DBs <- DBs %>% mutate(TrailMakingANNr_du_2_log=log(TrailMakingANNr_du_2+1))
describe(DBs$TrailMakingANNr_du_2_log)           #min 0.26 max 3.02
#Reverse log
DBs <- DBs %>% mutate(TrailMakingANNr_du_2_logr=3.02-TrailMakingANNr_du_2_log)
#Z score reverse log
DBs <- DBs %>% mutate(TrailMakingANNr_du_2_logrz = (TrailMakingANNr_du_2_logr - mean(TrailMakingANNr_du_2_logr, na.rm = TRUE))/sd(TrailMakingANNr_du_2_logr, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(TrailMakingANNr_du_2_st = ifelse(!is.na(TrailMakingANNr_du_2),1,0))

#Report
describe(DBs$TrailMakingANNr_du_2)                #min 0.3 max 19.45
describe(DBs$TrailMakingANNr_du_2_log)
describe(DBs$TrailMakingANNr_du_2_logr)
describe(DBs$TrailMakingANNr_du_2_logrz)
#table(DBs$TrailMakingANNr_du_2_st) 

#Figures
jpeg(filename = "64_TrailMakingANNr_du_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingANNr_du_2, breaks=20)         #skewed
dev.off() 
jpeg(filename = "64_TrailMakingANNr_du_2_log.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingANNr_du_2_log, breaks=20)
dev.off()
jpeg(filename = "64_TrailMakingANNr_du_2_logrz.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TrailMakingANNr_du_2_logrz, breaks=20)
dev.off() 


#7====SymbolDigit_2 
#Note: after excluding extremely high number of attempts  #-> 30899 (30925)
#Note: correction dropped because of severly skewed variable: Ratio of number of correct to number of attempts weighted by 12.5% (a better way to penalise attempts is needed, the number of attempts was not fixed, could press as many times as wanted)
table(DBs$SymbolDigit_c_2)
describe(DBs$SymbolDigit_c_2)                #min 0 max 36

#Z score 
DBs <- DBs %>% mutate(SymbolDigit_c_2_z = (SymbolDigit_c_2 - mean(SymbolDigit_c_2, na.rm = TRUE))/sd(SymbolDigit_c_2, na.rm = TRUE))
describe(DBs$SymbolDigit_c_2_z)
#Missing
#DBs <- DBs %>% mutate(SymbolDigit_c_2_st = ifelse(!is.na(SymbolDigit_c_2),1,0))

#Report  
describe(DBs$SymbolDigit_c_2)              #min 0 max 36
describe(DBs$SymbolDigit_c_2_z)
#table(DBs$SymbolDigit_2_st) 

#Figures 
jpeg(filename = "7_SymbolDigit_c_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$SymbolDigit_c_2, breaks=20)         #
dev.off() 
jpeg(filename = "7_SymbolDigit_c_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$SymbolDigit_c_2_z, breaks=20)         #
dev.off() 


#8====MatrixPattern_2 
#Note: Ratio of number correct over number viewed (f6374, the number of puzzles was fixed, so this reflects a true accuracy rate)
table(DBs$MatrixPattern_c_2)
table(DBs$MatrixPattern_2)
#Z score 
DBs <- DBs %>% mutate(MatrixPattern_2_z = (MatrixPattern_2 - mean(MatrixPattern_2, na.rm = TRUE))/sd(MatrixPattern_2, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(MatrixPattern_2_st = ifelse(!is.na(MatrixPattern_2),1,0))

#Report
describe(DBs$MatrixPattern_c_2)                #min 0 max 15
describe(DBs$MatrixPattern_2)                  #min 0 max 1
describe(DBs$MatrixPattern_2_z)
#table(DBs$MatrixPattern_2_st) 

#Figures
jpeg(filename = "8_MatrixPattern_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$MatrixPattern_2, breaks=20)            #normal
dev.off() 
jpeg(filename = "8_MatrixPattern_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$MatrixPattern_2_z, breaks=20)          #normal
dev.off() 


#9====TowerRearranging_2 
table(DBs$TowerRearranging_c_2)

#Z score 
DBs <- DBs %>% mutate(TowerRearranging_c_2_z = (TowerRearranging_c_2 - mean(TowerRearranging_c_2, na.rm = TRUE))/sd(TowerRearranging_c_2, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(TowerRearranging_c_2_st = ifelse(!is.na(TowerRearranging_c_2),1,0))

#Report
describe(DBs$TowerRearranging_c_2)                 #min 0 max 18
describe(DBs$TowerRearranging_c_2_z)
#table(DBs$TowerRearranging_c_2_st) 

#Figures
jpeg(filename = "9_TowerRearranging_c_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TowerRearranging_c_2, breaks=20)           #beautifully normal
dev.off() 
jpeg(filename = "9_TowerRearranging_c_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$TowerRearranging_c_2_z, breaks=20)         #normal (dentated)
dev.off() 


#10====PairedLearning_2
table(DBs$PairedLearning_2)

#Z score 
DBs <- DBs %>% mutate(PairedLearning_2_z = (PairedLearning_2 - mean(PairedLearning_2, na.rm = TRUE))/sd(PairedLearning_2, na.rm = TRUE))
#Missing
#DBs <- DBs %>% mutate(PairedLearning_2_st = ifelse(!is.na(PairedLearning_2),1,0))

#Report
describe(DBs$PairedLearning_2)                #min 0 max 10
describe(DBs$PairedLearning_2_z)
#table(DBs$PairedLearning_2_st) 

#Figures
jpeg(filename = "10_PairedLearning_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairedLearning_2, breaks=20)         #mildly negatively skewed
dev.off() 
jpeg(filename = "10_PairedLearning_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
hist(DBs$PairedLearning_2_z, breaks=20)         #
dev.off() 


#============================================g-Factor
##Read 
#CogGSet1a    5 components    @Instance-2      Executive Function score
DBg1 <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gfactor_and_updateddiag/CogGSet1a.rds")
names(DBg1)
describe(DBg1$f)
DBg1 <- rename(DBg1, gF5_2=f)
DBg1s=subset(DBg1, select= c(EID,gF5_2))
DBg1s <- rename(DBg1s, IID=EID)

#CogGSet1    5 components    @Instance-2      Imputed    Executive Function score
DBg1i <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gfactor_and_updateddiag/CogGSet1b.rds")
names(DBg1i)
describe(DBg1i$f)
DBg1i <- rename(DBg1i, gF5i_2=f)
DBg1is=subset(DBg1i, select= c(EID,gF5i_2))
DBg1is <- rename(DBg1is, IID=EID)

#CogGSet2a    9 components    @Instance-2      
DBg2 <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gfactor_and_updateddiag/CogGSet2a.rds")
names(DBg2)
describe(DBg2$f)
DBg2 <- rename(DBg2, gF9_2=f)
DBg2s=subset(DBg2, select= c(EID,gF9_2))
DBg2s <- rename(DBg2s, IID=EID)

#CogGSet2    9 components    @Instance-2      Imputed
DBg2i <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gfactor_and_updateddiag/CogGSet2b.rds")
names(DBg2i)
describe(DBg2i$f)
DBg2i <- rename(DBg2i, gF9i_2=f)
DBg2is=subset(DBg2i, select= c(EID,gF9i_2))
DBg2is <- rename(DBg2is, IID=EID)

#CogGSet3a    4 components    @Baseline      
DBg3 <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gfactor_and_updateddiag/CogGSet3a.rds")
names(DBg3)
describe(DBg3$f)
DBg3 <- rename(DBg3, gF4_0=f)
DBg3s=subset(DBg3, select= c(EID,gF4_0))
DBg3s <- rename(DBg3s, IID=EID)

##Merge
DBs <- left_join(DBs,DBg1s,by="IID")
DBs <- left_join(DBs,DBg1is,by="IID")
DBs <- left_join(DBs,DBg2s,by="IID")
DBs <- left_join(DBs,DBg2is,by="IID")
DBs <- left_join(DBs,DBg3s,by="IID")
str(DBs)
describe(DBs$gF5_2)   #36043 
describe(DBs$gF5i_2)  #50746
describe(DBs$gF9_2)   #35601
describe(DBs$gF9i_2)  #50746
describe(DBs$gF4_0)   #48741
#Inspecting distribution
hist(DBs$gF5i_2)
hist(DBs$gF5_2)                      
hist(DBs$gF9_2)                     
hist(DBs$gF9i_2)
hist(DBs$gF4_0)
#Remove unwanted data frames
rm(DBg1, DBg1s, DBg1i, DBg1is, DBg2, DBg2s, DBg2i, DBg2is, DBg3, DBg3s) 




#=======================================================================Exclusions
#====Missing structure: Original
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2,
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))


#====Ethnicity from CoreDataClin1 to exclude non-Whites
DB_Cov2 <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_pheno_ethn.rds")
table(DB_Cov2$Ethn_b)                              #0 472658, 1 28905
map(DB_Cov2, ~sum(is.na(.)))
#merge ethnicity 
dim(DBs)                                        #502411     66
DBs <- left_join(DBs,DB_Cov2,by="IID")
names(DBs)
dim(DBs)                                        #502411     68
DBs %>% count(Ethn_b)                           #0 472611, 1 28899, NA 901
#exclude non-Whites
DBs=subset(DBs, Ethn_b==0 & !is.na(Ethn_b)) 
DBs %>% count(Ethn_b)                           #0 472611
DBs=subset(DBs, select= -c(Ethn,Ethn_b))
names(DBs)

#Missing structure: non-Whites
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2, 
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))


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
table(DBs$exclu_Kinship)                              #n=442563  30048
#DBs<-filter(DBs, exclu_Kinship==0)
DBs=subset(DBs, exclu_Kinship==0)
table(DBs$exclu_Kinship)                              #n=442563
DBs=subset(DBs, select= -c(exclu_Kinship))
names(DBs)

#Missing structure: non-1st relatives
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2,
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))


#=====Poor heterozygosity/missingness
#Read ukb_sqc_v2.txt file
#DB_hm <- read.table(gzfile("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_sqc_v2.txt.gz"), header=TRUE)
DB_hm <- read.table(gzfile("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_sqc_v2.txt.gz"))
str(DB_hm)
names(DB_hm)
DB_hms <- DB_hm[c(19,25)] 
str(DB_hms)
table(DB_hms$V19)       #968      het.missing.outliers (poor-quality genotypes for these samples)
table(DB_hms$V25)       #407,219  used.in.pca.calculation
DB_hms <- rename(DB_hms, het.missing.outliers=V19)
#Read the .fam file to get the IID (the order of participants in the ukb_sqc_v2.txt file matches the order in the .fam file)
DB_id3 <- read.table("/data/Wolfson-PNU-dementia/genotype_hardcalls/ukb59138_cal_chr1_v2_s488264.fam")
str(DB_id3)
describe(DB_id3$V1)
DB_id <- DB_id3[c(1)] 
DB_id <- rename(DB_id, IID=V1)
str(DB_id)
#Past ukb_sqc_v2.txt and IID
rm(DB_hm, DB_id3)
dim(DB_id)                                         #488377   1
dim(DB_hms)                                        #488377   2
DB_hmid <- merge(DB_id,DB_hms)
DB_hmid <- cbind(DB_id,DB_hms)
names(DB_hmid)
dim(DB_hmid)                                        #488377      2
table(DB_hmid$het.missing.outliers)                 #968
#Saved as: ukb_SQC_outliers.rds
saveRDS(DB_hmid, file="/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_SQC_outliers.rds")
DB_hmid <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_SQC_outliers.rds")
DB_hmid=subset(DB_hmid, select= -c(V25))

#Merge het.missing.outliers to the outcome dataset
dim(DBs)                                        #442563     66
#.DBs <- merge(DBs,DB_hmid,by="IID")
DBs <- left_join(DBs,DB_hmid,by="IID")
names(DBs)
dim(DBs)                                        #442563     67
DBs %>% count(het.missing.outliers)                           #0 429097, 1 907, NA 12559
#exclude outliers
DBs=subset(DBs, het.missing.outliers==0 | is.na(het.missing.outliers)) 
DBs %>% count(het.missing.outliers)                           #0 429097, NA 12559
DBs=subset(DBs, select= -c(het.missing.outliers))
names(DBs)

#Missing structure: No poor heterozygosity/missingness
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2,
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))



#=====Covariates: PC 
##Read Cov_LF file
DBc <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_Cov_LF.rds")
DBcs=subset(DBc, select= c(IID,age,sex,centre,pc1,pc5,FEV1_b0,FVC_b0))
map(DBcs, ~sum(!is.na(.)))

DBcs <- DBcs %>% mutate(PCmis = ifelse(!is.na(pc1) & !is.na(pc5) & !is.na(age),0,1))
table(DBcs$PCmis)                               #488221  14241
DBcs <- DBcs %>% mutate(LFmis = ifelse(!is.na(FVC_b0),0,1))
table(DBcs$LFmis)                               #353285 149177
DBcs=subset(DBcs, select= c(IID,PCmis,LFmis))
str(DBcs)

#Merge Covariates&LF missing to the outcome dataset
dim(DBs)                                        #441656     66
DBs <- left_join(DBs,DBcs,by="IID")
names(DBs)
dim(DBs)                                        #441656     68
#exclude those without genetic data
DBs %>% count(PCmis)                            #0 429097, 1 12559
DBs=subset(DBs, PCmis==0) 
DBs %>% count(PCmis)                            #0 429097
DBs=subset(DBs, select= -c(PCmis))
names(DBs)

#Missing structure: No PC missing
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2,
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))


#=====Prevalent dementia
DBs <- rename(DBs, pDementia_0=ML_C42C240Xf41270_DementiaDiagyrBin0)
DBs <- rename(DBs, DementiaAge=ML_C42C240Xf41270_DementiaAge)
hist(DBs$DementiaAge)
describe(DBs$DementiaAge)           #min 38.9 max 84.15 n 6251
DBs <- rename(DBs, Age_2=R_Age_f21003_2)
hist(DBs$Age_2)
describe(DBs$Age_2)                 #min 45 max 83  n 47641         
DBs <- DBs %>% mutate(pDementia_2 = ifelse(DementiaAge<=Age_2,1,0))
table(DBs$pDementia_0)              #429006     91
table(DBs$pDementia_2)              #84 21

#exclude those with dementia at baseline               
DBs %>% count(pDementia_0)                           #429006     91
DBs=subset(DBs, pDementia_0==0) 
DBs %>% count(pDementia_0)                           #429006
DBs=subset(DBs, select= -c(pDementia_0))
names(DBs)

#exclude those with dementia at instance 2
DBs$PairsMatching_2[DBs$pDementia_2==1] <- NA
DBs$PairsMatching_2_logrz[DBs$pDementia_2==1] <- NA

DBs$ReactionTime_2[DBs$pDementia_2==1] <- NA
DBs$ReactionTime_2_logrz[DBs$pDementia_2==1] <- NA

DBs$ProspectiveMemory_2[DBs$pDementia_2==1] <- NA
DBs$ProspectiveMemory_2_br[DBs$pDementia_2==1] <- NA

DBs$FluidIntelligence_2[DBs$pDementia_2==1] <- NA
DBs$FluidIntelligence_2_z[DBs$pDementia_2==1] <- NA

DBs$NumericMemory_2[DBs$pDementia_2==1] <- NA
DBs$NumericMemory_2_z[DBs$pDementia_2==1] <- NA

DBs$TrailMakingN_du_2[DBs$pDementia_2==1] <- NA
DBs$TrailMakingN_du_2_logrz[DBs$pDementia_2==1] <- NA

DBs$TrailMakingAN_du_2[DBs$pDementia_2==1] <- NA
DBs$TrailMakingAN_du_2_logrz[DBs$pDementia_2==1] <- NA

DBs$TrailMakingANN_du_2[DBs$pDementia_2==1] <- NA
DBs$TrailMakingANN_du_2_rz[DBs$pDementia_2==1] <- NA

DBs$SymbolDigit_c_2[DBs$pDementia_2==1] <- NA
DBs$SymbolDigit_c_2_z[DBs$pDementia_2==1] <- NA

DBs$MatrixPattern_2[DBs$pDementia_2==1] <- NA
DBs$MatrixPattern_2_z[DBs$pDementia_2==1] <- NA

DBs$TowerRearranging_c_2[DBs$pDementia_2==1] <- NA
DBs$TowerRearranging_c_2_z[DBs$pDementia_2==1] <- NA

DBs$PairedLearning_2[DBs$pDementia_2==1] <- NA
DBs$PairedLearning_2_z[DBs$pDementia_2==1] <- NA

DBs$gF4_0[DBs$pDementia_2==1] <- NA
DBs$gF5_2[DBs$pDementia_2==1] <- NA
DBs$gF5i_2[DBs$pDementia_2==1] <- NA
DBs$gF9_2[DBs$pDementia_2==1] <- NA
DBs$gF9i_2[DBs$pDementia_2==1] <- NA

#Missing structure: No dementia before assessments
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0_logrz, PairsMatching_2_logrz, 
                             ReactionTime_0_logrz, ReactionTime_2_logrz, 
                             ProspectiveMemory_0_br, ProspectiveMemory_2_br, 
                             FluidIntelligence_0_z, FluidIntelligence_2_z, 
                             NumericMemory_0_z, NumericMemory_2_z,
                             TrailMakingN_du_2_logrz, TrailMakingAN_du_2_logrz, 
                             TrailMakingANN_du_2_rz, SymbolDigit_c_2_z, 
                             MatrixPattern_2_z, TowerRearranging_c_2_z, 
                             PairedLearning_2_z, gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))


#=====Saved as: ukb_pheno_cogTS.rds
print("Saved as: ukb_pheno_cogTS.rds")
saveRDS(DBs, file="/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_pheno_cogTS.rds")
print("Saved as: ukb_pheno_cogTS.txt")
write.table(DBs, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_pheno_cogTS.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#=====Covariates: LF
#exclude those without LF data              --->>> NOTE: Decided NOT to exclude 
DBs %>% count(LFmis)                           #0 322887, 1 106119
DBs=subset(DBs, LFmis==0) 
DBs %>% count(LFmis)                           #0 322887
DBs=subset(DBs, select= -c(LFmis))
names(DBs)
DBss <- subset(DBs,select= c(FID, IID, PairsMatching_0, PairsMatching_2, ReactionTime_0, 
                             ReactionTime_2, ProspectiveMemory_0, ProspectiveMemory_2, 
                             FluidIntelligence_0, FluidIntelligence_2, NumericMemory_0, 
                             NumericMemory_2, TrailMakingN_du_2, TrailMakingAN_du_2, 
                             TrailMakingANN_du_2, SymbolDigit_c_2, MatrixPattern_2, 
                             TowerRearranging_c_2, PairedLearning_2,
                             gF4_0, gF5_2, gF5i_2, gF9_2, gF9i_2))
map(DBss, ~sum(!is.na(.)))

#Remove unwanted data frames
rm(DB_Cov2, DB_id, DB_hmid, DB_hms, DBc, DBcs, Exclu_kinship, kinship)



#=====Exploring correlations (Graphs in Corr_Cog file)
# @Baseline
DBcor <- subset(DBs,select= c(PairsMatching_0_log, ReactionTime_0_log, FluidIntelligence_0,  
                              NumericMemory_0, gF4_0))
DBcor <- rename(DBcor, 
                Pairs_Matching = 'PairsMatching_0_log',
                Reaction_Time = 'ReactionTime_0_log',
                Fluid_Intelligence = 'FluidIntelligence_0',
                Numeric_Memory = 'NumericMemory_0' ,
                gFactor_4items = 'gF4_0')
Cor <- cor(DBcor, method = "pearson", use = "complete.obs")
Cor <- round(Cor, 2)

#Figures
library(corrplot)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1, number.cex=1)


# @Instance 2
DBcor <- subset(DBs,select= c(PairsMatching_2_log, ReactionTime_2_log, FluidIntelligence_2,  
                              NumericMemory_2, TrailMakingANN_du_2, SymbolDigit_c_2,  
                              MatrixPattern_2, TowerRearranging_c_2, PairedLearning_2,
                              gF5i_2, gF9i_2))
DBcor <- rename(DBcor, 
                Pairs_Matching = 'PairsMatching_2_log',
                Reaction_Time = 'ReactionTime_2_log',
                Fluid_Intelligence = 'FluidIntelligence_2',
                Numeric_Memory = 'NumericMemory_2' ,
                Trail_Making_BA = 'TrailMakingANN_du_2',
                Symbol_Digit = 'SymbolDigit_c_2',
                Matrix_Pattern = 'MatrixPattern_2',
                Tower_Rearranging = 'TowerRearranging_c_2',
                Paired_Learning = 'PairedLearning_2',
                gFactor_5items = 'gF5i_2',
                gFactor_9items = 'gF9i_2')
Cor <- cor(DBcor, method = "pearson", use = "complete.obs")
Cor <- round(Cor, 2)

#Figures
library(corrplot)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9, number.cex=0.8)

rm(Cor,DBcor)



###============================================Outcome files
##Read extracted file, managed
DBs <- readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ukb_pheno_cogTS.rds")

#1.Pairs matching file
## @B
PairsMatching = DBs%>% select(contains("FID"), contains("IID"), PairsMatching_0_logrz)
map(PairsMatching, ~sum(!is.na(.)))
PairsMatching <- PairsMatching[complete.cases(PairsMatching),]
print("Structure- Reaction Time at instance 0:")
str(PairsMatching)
write.table(PairsMatching, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/PairsMatching_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: PairsMatching_0.txt")
## @I2
PairsMatching = DBs%>% select(contains("FID"), contains("IID"), PairsMatching_2_logrz)
map(PairsMatching, ~sum(!is.na(.)))
PairsMatching <- PairsMatching[complete.cases(PairsMatching),]
print("Structure- Reaction Time at instance 2:")
str(PairsMatching)
write.table(PairsMatching, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/PairsMatching_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: PairsMatching_2.txt")

#2.Reaction Time file
## @B
ReactionTime = DBs%>% select(contains("FID"), contains("IID"), ReactionTime_0_logrz)
map(ReactionTime, ~sum(!is.na(.)))
ReactionTime <- ReactionTime[complete.cases(ReactionTime),]
print("Structure- Reaction Time at instance 0:")
str(ReactionTime)
saveRDS(ReactionTime, file="/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ReactionTime_0.rds")
print("Saved as: ReactionTime_0.rds")
write.table(ReactionTime, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ReactionTime_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: ReactionTime_0.txt")
## @I2
ReactionTime = DBs%>% select(contains("FID"), contains("IID"), ReactionTime_2_logrz)
map(ReactionTime, ~sum(!is.na(.)))
ReactionTime <- ReactionTime[complete.cases(ReactionTime),]
print("Structure- Reaction Time at instance 2:")
str(ReactionTime)
saveRDS(ReactionTime, file="/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ReactionTime_2.rds")
print("Saved as: ReactionTime_2.rds")
write.table(ReactionTime, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/ReactionTime_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: ReactionTime_2.txt")

#4.Fluid Intelligence file
## @B
FluidIntelligence = DBs%>% select(contains("FID"), contains("IID"), FluidIntelligence_0_z)
map(FluidIntelligence, ~sum(!is.na(.)))
FluidIntelligence <- FluidIntelligence[complete.cases(FluidIntelligence),]
print("Structure- Reaction Time at instance 0:")
str(FluidIntelligence)
write.table(FluidIntelligence, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/FluidIntelligence_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: FluidIntelligence_0.txt")
## @I2
FluidIntelligence = DBs%>% select(contains("FID"), contains("IID"), FluidIntelligence_2_z)
map(FluidIntelligence, ~sum(!is.na(.)))
FluidIntelligence <- FluidIntelligence[complete.cases(FluidIntelligence),]
print("Structure- Reaction Time at instance 2:")
str(FluidIntelligence)
write.table(FluidIntelligence, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/FluidIntelligence_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: FluidIntelligence_2.txt")

#5.Numeric Memory file
## @B
NumericMemory = DBs%>% select(contains("FID"), contains("IID"), NumericMemory_0_z)
map(NumericMemory, ~sum(!is.na(.)))
NumericMemory <- NumericMemory[complete.cases(NumericMemory),]
print("Structure- Reaction Time at instance 0:")
str(NumericMemory)
write.table(NumericMemory, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/NumericMemory_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: NumericMemory_0.txt")
## @I2
NumericMemory = DBs%>% select(contains("FID"), contains("IID"), NumericMemory_2_z)
map(NumericMemory, ~sum(!is.na(.)))
NumericMemory <- NumericMemory[complete.cases(NumericMemory),]
print("Structure- Reaction Time at instance 2:")
str(NumericMemory)
write.table(NumericMemory, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/NumericMemory_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: NumericMemory_2.txt")

#6.TMT B-A file
## @I2
TrailMaking = DBs%>% select(contains("FID"), contains("IID"), TrailMakingANN_du_2_rz)
map(TrailMaking, ~sum(!is.na(.)))
TrailMaking <- TrailMaking[complete.cases(TrailMaking),]
print("Structure- Reaction Time at instance 2:")
str(TrailMaking)
write.table(TrailMaking, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/TrailMaking_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: TrailMaking_2.txt")

#7.Symbol Digit file
## @I2
SymbolDigit = DBs%>% select(contains("FID"), contains("IID"), SymbolDigit_c_2_z)
map(SymbolDigit, ~sum(!is.na(.)))
SymbolDigit <- SymbolDigit[complete.cases(SymbolDigit),]
print("Structure- Reaction Time at instance 2:")
str(SymbolDigit)
write.table(SymbolDigit, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/SymbolDigit_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: SymbolDigit_2.txt")

#8.Matrix Pattern file
## @I2
MatrixPattern = DBs%>% select(contains("FID"), contains("IID"), MatrixPattern_2_z)
map(MatrixPattern, ~sum(!is.na(.)))
MatrixPattern <- MatrixPattern[complete.cases(MatrixPattern),]
print("Structure- Reaction Time at instance 2:")
str(MatrixPattern)
write.table(MatrixPattern, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/MatrixPattern_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: MatrixPattern_2.txt")

#9.Tower Rearranging file
## @I2
TowerRearranging = DBs%>% select(contains("FID"), contains("IID"), TowerRearranging_c_2_z)
map(TowerRearranging, ~sum(!is.na(.)))
TowerRearranging <- TowerRearranging[complete.cases(TowerRearranging),]
print("Structure- Reaction Time at instance 2:")
str(TowerRearranging)
write.table(TowerRearranging, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/TowerRearranging_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: TowerRearranging_2.txt")

#10.Paired Learning file
## @I2
PairedLearning = DBs%>% select(contains("FID"), contains("IID"), PairedLearning_2_z)
map(PairedLearning, ~sum(!is.na(.)))
PairedLearning <- PairedLearning[complete.cases(PairedLearning),]
print("Structure- Reaction Time at instance 2:")
str(PairedLearning)
write.table(PairedLearning, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/PairedLearning_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: PairedLearning")

#11.g-Factor @B
gF4_0 = DBs%>% select(contains("FID"), contains("IID"), gF4_0)
map(gF4_0, ~sum(!is.na(.)))
gF4_0 <- gF4_0[complete.cases(gF4_0),]
str(gF4_0)
write.table(gF4_0, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gF4_0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#12.g-Factor Exec Func @I2
gF5i_2 = DBs%>% select(contains("FID"), contains("IID"), gF5i_2)
map(gF5i_2, ~sum(!is.na(.)))
gF5i_2 <- gF5i_2[complete.cases(gF5i_2),]
str(gF5i_2)
write.table(gF5i_2, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gF5i_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#13.g-Factor Full @I2
gF9i_2 = DBs%>% select(contains("FID"), contains("IID"), gF9i_2)
map(gF9i_2, ~sum(!is.na(.)))
gF9i_2 <- gF9i_2[complete.cases(gF9i_2),]
str(gF9i_2)
write.table(gF9i_2, file = "/data/Wolfson-PNU-dementia/UKB_Projects/lungdev/gF9i_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

