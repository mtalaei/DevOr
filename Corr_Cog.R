#-------------------------------------------------------------------------------
# Correlation between Cognitive Function measures (instance 0 @Baseline AND @2)
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) 

##====Read
DBs <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_cogTS.rds")

##====Set working directory for outputs (figures)
setwd("/data/home/hmy431/outputs/Figures/")


#=====Correlations before reversing and Z score
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
DBcor <- round(Cor, 2)

#Figures
library(corrplot)
jpeg(filename = "CF-GF_Corr_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1)
dev.off() 

jpeg(filename = "CF-GF_Corr_r_0.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1, number.cex=1)
dev.off() 

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
DBcor <- round(Cor, 2)

#Figures
library(corrplot)
jpeg(filename = "CF-GF_Corr_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9)
dev.off() 

jpeg(filename = "CF-GF_Corr_r_2.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9, number.cex=0.8)
dev.off() 

rm(Cor,DBcor)


#=====Correlations for CF items at their final form 

# @Baseline
#PairsMatching_0_logrz, ReactionTime_0_logrz, FluidIntelligence_0_z, NumericMemory_0_z

DBcor <- subset(DBs,select= c(PairsMatching_0_logrz, ReactionTime_0_logrz, FluidIntelligence_0_z, NumericMemory_0_z, gF4_0))
DBcor <- rename(DBcor, 
                Pairs_Matching = 'PairsMatching_0_logrz',
                Reaction_Time = 'ReactionTime_0_logrz',
                Fluid_Intelligence = 'FluidIntelligence_0_z',
                Numeric_Memory = 'NumericMemory_0_z',
                gFactor_4items = 'gF4_0')
Cor <- cor(DBcor, method = "pearson", use = "complete.obs")
Cor <- round(Cor, 2)

#Figures
library(corrplot)
jpeg(filename = "CF-GF_Corr_0_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1)
dev.off() 

jpeg(filename = "CF-GF_Corr_r_0_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 1, number.cex=1)
dev.off() 

# @Instance 2
#TrailMakingANN_du_2_rz, SymbolDigit_c_2_z, MatrixPattern_2_z, TowerRearranging_c_2_z, PairedLearning_2_z

#DBcor <- subset(DBs,select= c(TrailMakingANN_du_2_rz, SymbolDigit_c_2_z, MatrixPattern_2_z,
#                              TowerRearranging_c_2_z, PairedLearning_2_z, gF5i_2, gF9i_2))
DBcor <- subset(DBs,select= c(PairsMatching_2_logrz, ReactionTime_2_logrz, FluidIntelligence_2_z, 
                              NumericMemory_2_z, TrailMakingANN_du_2_rz, SymbolDigit_c_2_z, 
                              MatrixPattern_2_z, TowerRearranging_c_2_z, PairedLearning_2_z, gF5i_2, gF9i_2))
DBcor <- rename(DBcor, 
                Pairs_Matching = 'PairsMatching_2_logrz',
                Reaction_Time = 'ReactionTime_2_logrz',
                Fluid_Intelligence = 'FluidIntelligence_2_z',
                Numeric_Memory = 'NumericMemory_2_z' ,
                Trail_Making_BA = 'TrailMakingANN_du_2_rz',
                Symbol_Digit = 'SymbolDigit_c_2_z',
                Matrix_Pattern = 'MatrixPattern_2_z',
                Tower_Rearranging = 'TowerRearranging_c_2_z',
                Paired_Learning = 'PairedLearning_2_z',
                gFactor_5items = 'gF5i_2',
                gFactor_9items = 'gF9i_2')
Cor <- cor(DBcor, method = "pearson", use = "complete.obs")
Cor <- round(Cor, 2)

#Figures
library(corrplot)
jpeg(filename = "CF-GF_Corr_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9)
dev.off() 

jpeg(filename = "CF-GF_Corr_r_2_z.jpeg", width = 8, height = 8, units = 'in', res = 300)
corrplot(Cor, method = 'number', type = "upper", tl.col = "black", tl.srt = 45, tl.cex = 0.9, number.cex=0.8)
dev.off() 

rm(Cor,DBcor)