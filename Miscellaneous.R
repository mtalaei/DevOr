#===========================================================
#                 Miscellaneous   @HPC
#===========================================================
rm(list = ls(all.names = TRUE)) 

#----@HPC
library(dplyr)
library(tidyverse)
library(data.table)
library(R.utils)
library("tidyr")
#Set working directory
setwd("/data/home/hmy431/")

#List
#========Manhattan plot
#========Comparing with look-up
#========Ast-MI
#========BP variables
#========MAF min-max
#========Comparing effect alleles with Laura output
#========Replicating Laura sample size
#========DaMa Ethn
#========Extracting summary statistics from Regenie outputs 
  #Look-up for rs11865499 & rs1978487 in meta-GWAS for AD
  #for Pairs Matching
  #for LF
  #All (--> Coloclized_Variants)
#========Providing complete genetic file including position
#========Investigating mismatch in merging meta-GWAS with UKB
#========Extracting summary statistics from Regenie outputs




##============================================Manhattan plot
library(qqman)
DB <- read.table("outputs/FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DB$LOG10P <- as.numeric(DB$LOG10P)
DB$Pval <- 10^(DB$LOG10P*-1)
names(DB)
#Plot
manhattan(DB, chr="CHROM", bp="GENPOS", snp="ID", p="Pval" )
snpsOfInterest <- c("rs9267531", "rs237833", "rs548092276", "rs11275011", "rs28496034")
snpsOfInterest
manhattan(DB, chr="CHROM", bp="GENPOS", snp="ID", p="Pval", highlight = snpsOfInterest, annotatePval = 0.001)
qq(DB$Pval)



##============================================Comparing with look-up
DB <- read.table("outputs/FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DBs<-DB[match(c("rs11229063", "rs28446321", "rs665058", "rs972936"),DB$ID),]
#P-value
DBs$LOG10P <- as.numeric(DBs$LOG10P)
DBs$Pval <- 10^(DBs$LOG10P*-1)
DBs<-subset(DBs, select=c(ID, ALLELE0, ALLELE1, BETA, Pval))
DBs


##============================================Ast-MI
#Note: the old version volume is 187M (now 361M)
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/UKBW_52532_110822A_CogDem.rds")
DBs = DB%>% select(contains("EID"), 
                   contains("ML_C42C240Xf41270f20002_Asthma"), 
                   contains("ML_C42C240Xf41270f20002_IschaemicHeartMI"),
                   contains("ML_C42C240Xf41270f20002_Angina"),
                   contains("ML_f20004_OP1069Heart"),
                   contains("ML_f20004_OP1070Angioplasty"),
                   contains("ML_f20004_OP1095CoronaryArteryBypass"),
                   contains("R_ML_IschaemicHeartMIComposite_f6150f3894C42C240Xf41270f20002"))
DBss = DBs%>% select("EID", 
                    "ML_C42C240Xf41270f20002_Asthma",
                    "ML_C42C240Xf41270f20002_AsthmaDiagyrBin",
                    "ML_C42C240Xf41270f20002_IschaemicHeartMI",
                    "ML_C42C240Xf41270f20002_IschaemicHeartMIDiagyrBin",
                    "R_ML_IschaemicHeartMIComposite_f6150f3894C42C240Xf41270f20002_0")
DBss = DBss%>% rename(IID=EID)
DBss = DBss%>% rename(Ast=ML_C42C240Xf41270f20002_Asthma)
DBss = DBss%>% rename(AstBin=ML_C42C240Xf41270f20002_AsthmaDiagyrBin)
DBss = DBss%>% rename(MI=ML_C42C240Xf41270f20002_IschaemicHeartMI)
DBss = DBss%>% rename(MIBin=ML_C42C240Xf41270f20002_IschaemicHeartMIDiagyrBin)
DBss = DBss%>% rename(IHD_0=R_ML_IschaemicHeartMIComposite_f6150f3894C42C240Xf41270f20002_0)

#Save
export(DBss,"/data/Wolfson-PNU-dementia/lungdev/ukb_AstMI.dta")


#============================================BP variables
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/UKBW_52532_110822A_CogDem.rds")
##Select components
DB_H = DB%>% select("EID", 
                    "R_Systolic_f4080_mean_0",
                    "systolic_blood_pressure_automated_reading_f4080_0_0",
                    "systolic_blood_pressure_automated_reading_f4080_0_1",
                    "R_Diastolic_f4079_mean_0",
                    "diastolic_blood_pressure_automated_reading_f4079_0_0",
                    "diastolic_blood_pressure_automated_reading_f4079_0_1",
                    "R_MedBloodPressureComposite_f6153f6177_0",
                    "ML_f41270_I10Hypertension", 
                    "ML_f20002_SR1072EssentialHypertension", 
                    "ML_C2409_I10essentialprimaryHypertension",
                    "ML_C2409_I10essentialprimaryHypertensionAge",
                    "R_Age_f21003_0")
names(DB_H)
#Defining HTN at baseline
#Note: mutate(R_Systolic_f4080_mean_0=rowMeans(cbind(systolic_blood_pressure_automated_reading_f4080_0_0,systolic_blood_pressure_automated_reading_f4080_0_1), na.rm = TRUE))

describe(DB_H$R_Systolic_f4080_mean_0)                                    #min 65 max 268  n 472294
describe(DB_H$R_Diastolic_f4079_mean_0)                                   #min 32 max 147.5  n 472299

describe(DB_H$systolic_blood_pressure_automated_reading_f4080_0_1)        #min 56 max 259  n 461197
describe(DB_H$diastolic_blood_pressure_automated_reading_f4079_0_1)       #min 30 max 147  n 461201

#Rename
DB_H = DB_H%>% rename(SBPm_0=R_Systolic_f4080_mean_0)
DB_H = DB_H%>% rename(DBPm_0=R_Diastolic_f4079_mean_0)
DB_H = DB_H%>% rename(IID=EID)

#Select
DB_H = DB_H%>% select(IID, SBPm_0, DBPm_0)

#Save
export(DB_H,"/data/Wolfson-PNU-dementia/lungdev/ukb_BP0.dta")
                    

#============================================MAF min-max
DB <- read.table("/data/home/hmy431/data/Variants.txt", header=TRUE)
summary(DB$MAF)
hist(DB$MAF)
#5%
DB$MAF05 <- cut(DB$MAF, breaks=c(0.000058, 0.04999, 0.5), labels=c('<5%', '>=5%'))
table(DB$MAF05)               #<5% 4,798    >5% 10,842
sum(table(DB$MAF05))          #15,640
TT <-table(DB$MAF05)
prop.table(TT)                #0.3067775 0.6932225
#1%
DB$MAF01 <- cut(DB$MAF, breaks=c(0.000058, 0.00999, 0.5), labels=c('<1%', '>=1%'))
table(DB$MAF01)               #<5% 4,798    >5% 10,842
TT <-table(DB$MAF01)
prop.table(TT)                #0.001023018 0.998976982

#KAT8
genelist <- DB %>% count(GENE)
KAT8 <- subset(DB, GENE=="KAT8")
summary(KAT8$MAF)
table(KAT8$MAF05) 

#============================================Comparing effect alleles with Laura output
#Laura's new file containing alleles
DB <- read.table("/data/home/hmy431/outputs/fromLaura/LF_sumstats.txt", header=TRUE)
head(DB)
DB <- subset(DB, select=c(ID, ALLELE0, ALLELE1, beta_fvc, beta_ratio))
DB <- rename(DB, "ALLELE0_LF"="ALLELE0")
DB <- rename(DB, "ALLELE1_LF"="ALLELE1")
str(DB)

#variants with the signal in colocalization 
DBr <- readRDS("/data/home/hmy431/data/Coloclized_Variants.rds")
head(DBr)
DBr <- subset(DBr, select=c(ID, ALLELE0, ALLELE1, BETA, CF))
DBr <- rename(DBr, "ALLELE0_CF"="ALLELE0")
DBr <- rename(DBr, "ALLELE1_CF"="ALLELE1")
DBr <- rename(DBr, "BETA_CF"="BETA")
str(DBr)

#merging
DBs <- left_join(DBr,DB,by="ID")
DBs


#============================================Replicating Laura sample size
DB <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")      #raw extracted file
DBs <- subset(DB, select=c(FID, IID, FEV1_b0, FVC_b0, FEV1FVC_b0))
names(DBs)
map(DBs, ~sum(!is.na(.)))         #502462
#PC
DBc <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_Cov_LF.rds")
DBcs=subset(DBc, select= c(IID,pc1,pc5))
DBcs <- DBcs %>% mutate(PCmis = ifelse(!is.na(pc1) & !is.na(pc5),0,1))
table(DBcs$PCmis)                               #488221  14241
DBcss=subset(DBcs, select= c(IID,PCmis))
DBs <- left_join(DBs,DBcss,by="IID")
names(DBs)
DBs %>% count(PCmis)                            #0 488221, 1 14241
DBs=subset(DBs, PCmis==0) 
DBs %>% count(PCmis)                            #0 488221
DBs=subset(DBs, select= -c(PCmis))
map(DBs, ~sum(!is.na(.)))         #488221
#White British 
DB_Cov2 <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_ethn.rds")
DB_Cov2 <- DB_Cov2 %>% mutate(BWhite = ifelse(Ethn==1001, 0, 1))
table(DB_Cov2$BWhite)                            #0 442552  1 59011
DBs <- left_join(DBs,DB_Cov2,by="IID")
DBs=subset(DBs, BWhite==0 & !is.na(BWhite)) 
DBs %>% count(BWhite)                           #0 430985
map(DBs, ~sum(!is.na(.)))
#Poor heterozygosity/missingness
DB_hmid <- readRDS("/data/Wolfson-PNU-dementia/lungdev/ukb_SQC_outliers.rds") #Read ukb_SQC_outliers.rds
table(DB_hmid$het.missing.outliers)                #968
DB_hmid=subset(DB_hmid, select= -c(V25))
DBs <- left_join(DBs,DB_hmid,by="IID")
names(DBs)
DBs %>% count(het.missing.outliers)                           #0 430182, 1 803
DBs=subset(DBs, het.missing.outliers==0) 
DBs %>% count(het.missing.outliers)                           #0 430182
map(DBs, ~sum(!is.na(.)))
#Kinship (Relatedness)
kinship = read_table("/data/Wolfson-PNU-dementia/UKB/ukb_helper_files/ukb59138_rel_s488264.dat")
kinship$Kinship_0.125 <- ifelse(kinship$Kinship >= 0.125, c("higher"), c("lower"))  #4th degree
table(kinship$Kinship_0.125)                            #n=73090/34026
kinship$Kinship_0.25 <- ifelse(kinship$Kinship >= 0.25, c("higher"), c("lower"))    #3rd degree
table(kinship$Kinship_0.25)
Exclu_kinship = kinship %>% filter(Kinship>0.125) %>% select(ID1) %>% rename("IID"="ID1") 
DBs <- DBs  %>% mutate(exclu_Kinship = ifelse(IID %in% Exclu_kinship$IID, "1", "0"))
table(DBs$exclu_Kinship)                              #0 418950  1 11232
DBs=subset(DBs, exclu_Kinship==0)
table(DBs$exclu_Kinship)                              #n=418950
map(DBs, ~sum(!is.na(.)))



#============================================DaMa Ethn
#Ethnicity from CoreDataClin1 to exclude non-Whites (codes found in zArchive_DaMa_cog_2.R)
#DB <- read_tsv("/data/Wolfson-PNU-dementia/datasets/ukb_extracted/47532_Coded/1_ukbpheno_CoreDataClin1.tsv")
DB <- read_tsv("/data/Wolfson-PNU-dementia/UKB/datasets/ukb_extracted/59138r47532_Coded/1_ukbpheno_CoreDataClin1.tsv")
DB_Cov2 = DB%>% select(contains("EID"),
                       contains("Ethnic background.0.0"))
DB_Cov2 <- rename(DB_Cov2, "IID"="EID")
DB_Cov2 <- rename(DB_Cov2, "Ethn"="Ethnic background.0.0")
str(DB_Cov2)
print("----Ethnic background")
table(DB_Cov2$Ethn)
map(DB_Cov2, ~sum(is.na(.)))
DB_Cov2 <- DB_Cov2 %>% mutate(Ethn_b = ifelse(Ethn==1 | Ethn==1001 | Ethn==1002 | Ethn==1003, 0, 1))
saveRDS(DB_Cov2, file="/data/Wolfson-PNU-dementia/lungdev/ukb_pheno_ethn.rds")

#merge ethnicity   
DB_CF <- merge(DB_CF,DB_Cov2,by="IID")
str(DB_CF)
dim(DB_CF)                                        #471340
DB_CF %>% count(Ethn_b)                           #0 442606, 1 27862, NA 872
#exclude non-Whites
DB_CF=subset(DB_CF, Ethn_b==0 & !is.na(Ethn_b)) 
DB_CF %>% count(Ethn_b)                           #0 442606
DB_CF=subset(DB_CF, select= -c(Ethn,Ethn_b))
map(DB_CF, ~sum(!is.na(.)))

#==White British (Exploratory) 
DB_Cov2 <- DB_Cov2 %>% mutate(BWhite = ifelse(Ethn==1001, 0, 1))
table(DB_Cov2$BWhite)
DBs %>% count(BWhite)
table(DBs$BWhite, DBs$Ethn_b)
#exclude non-Whites
DBs=subset(DBs, BWhite==0 & !is.na(BWhite)) 
DBs %>% count(BWhite)                           #0 472658
DBs=subset(DBs, select= -c(Ethn,Ethn_b,BWhite))
names(DBs)
#Missing structure: non-Whites
map(DBs, ~sum(!is.na(.)))


#============================================Extracting summary statistics from Regenie outputs 
#========Look-up for rs11865499 & rs1978487 in meta-GWAS for AD
DB <- read.table("/data/home/hmy431/data/ukb_GWASAD_add.txt", header=TRUE)
DB1<-DB[match(c("rs11865499", "rs1978487"),DB$ID),]
DB1
DB_LF <- read.table("/data/home/hmy431/outputs/fromLaura/LF_sumstats.txt", header=TRUE)
DB2<-DB_LF[match(c("rs11865499", "rs1978487"),DB_LF$ID),]
DB2

#========for Pairs Matching
setwd("/data/home/hmy431/outputs/")  

DB <- read.table("PairsMatching_0_H/ukb_step2_Pairs_0.regenie", header=TRUE)
DB1<-DB[match(c("rs138259061"),DB$ID),]
DB1$P_Pairs <- 10^(DB1$LOG10P*-1)
str(DB1)


#========for LF
setwd("/data/home/hmy431/outputs/")  

#FVC
DB <- read.table("FVC_0/ukb_step2_FVC_0.regenie", header=TRUE)
DB1<-DB[match(c("rs9267531","rs237833","rs548092276","rs11275011","rs2084448",
                "rs138259061","rs1978487","rs11865499","rs113154802","rs539078574",
                "rs28496034","rs2297086","rs75614054","rs6120880"),DB$ID),]
str(DB1)
#P-value
DB1$LOG10P <- as.numeric(DB1$LOG10P)
DB1$P_FVC <- 10^(DB1$LOG10P*-1)
#Select
DB1=subset(DB1, select= c(ID, ALLELE0, ALLELE1, BETA, P_FVC))
DB1 <- rename(DB1, ALLELE0_FVC=ALLELE0)
DB1 <- rename(DB1, ALLELE1_FVC=ALLELE1)
DB1 <- rename(DB1, BETA_FVC=BETA)

#Ratio
DB <- read.table("Ratio_0/ukb_step2_Ratio_0.regenie", header=TRUE)
DB2<-DB[match(c("rs9267531","rs237833","rs548092276","rs11275011","rs2084448",
                "rs138259061","rs1978487","rs11865499","rs113154802","rs539078574",
                "rs28496034","rs2297086","rs75614054","rs6120880"),DB$ID),]
str(DB2)
#P-value
DB2$LOG10P <- as.numeric(DB2$LOG10P)
DB2$P_Ratio <- 10^(DB2$LOG10P*-1)
#Select
DB2=subset(DB2, select= c(ID, ALLELE0, ALLELE1, BETA, P_Ratio))
DB2 <- rename(DB2, ALLELE0_Ratio=ALLELE0)
DB2 <- rename(DB2, ALLELE1_Ratio=ALLELE1)
DB2 <- rename(DB2, BETA_Ratio=BETA)

#Merge FVC and Ratio
DB_LF <- merge(DB1,DB2,by="ID")
str(DB_LF)

#Reading CF data
DB <- readRDS("/data/home/hmy431/data/Coloclized_Variants.rds")
str(DB)
DB$LOG10P <- as.numeric(DB$LOG10P)
DB$P_CF <- 10^(DB$LOG10P*-1)
DB_CF=subset(DB, select= c(ID, ALLELE0, ALLELE1, BETA, P_CF, CF))
str(DB_CF)

#Merging CF with LF
DB <- left_join(DB_CF, DB_LF, by = "ID")
str(DB)
write_xlsx(DB,"/data/home/hmy431/data/CF_LF_Regenie.xlsx", col_names = TRUE, format_headers = TRUE)



#========All (--> Coloclized_Variants)
rm(list = ls(all.names = TRUE)) 
setwd("/data/home/hmy431/outputs/")     

#Fluid intelligence
DB <- read.table("FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DB_<-DB[match(c("rs9267531", "rs3117578", "rs575597758"),DB$ID),]
DB_<-DB_%>%mutate(CF="Fluid intelligence")
DBs <- DB_
#Pairs matching
DB <- read.table("PairsMatching_0/ukb_step2_Pairs_0.regenie", header=TRUE)
DB_<-DB[match(c("rs9267531", "rs3117578", "rs575597758"),DB$ID),]
DB_<-DB_%>%mutate(CF="Pairs matching")
DBs <- rbind(DBs, DB_)
#Fluid intelligence
DB <- read.table("FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DB_<-DB[match(c("rs237833", "rs548092276", "rs11275011"),DB$ID),]
DB_<-DB_%>%mutate(CF="Fluid intelligence")
DBs <- rbind(DBs, DB_)
#Pairs matching
DB <- read.table("PairsMatching_0/ukb_step2_Pairs_0.regenie", header=TRUE)
DB_<-DB[match(c("rs2084448", "rs138259061"),DB$ID),]
DB_<-DB_%>%mutate(CF="Pairs matching")
DBs <- rbind(DBs, DB_)

#Reaction time
DB <- read.table("ReactionTime_0/ukb_step2_RT_0.regenie", header=TRUE)
DB_<-DB[match(c("rs113154802", "rs539078574", "rs2227603", "rs2227592"),DB$ID),]
DB_<-DB_%>%mutate(CF="Reaction time")
DBs <- rbind(DBs, DB_)
#Fluid intelligence
DB <- read.table("FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DB_<-DB[match(c("rs113154802", "rs28496034"),DB$ID),]
DB_<-DB_%>%mutate(CF="Fluid intelligence")
DBs <- rbind(DBs, DB_)
#Reaction time
DB <- read.table("ReactionTime_0/ukb_step2_RT_0.regenie", header=TRUE)
DB_<-DB[match(c("rs2297086", "rs539078574", "rs75614054"),DB$ID),]
DB_<-DB_%>%mutate(CF="Reaction time")
DBs <- rbind(DBs, DB_)
#Fluid intelligence
DB <- read.table("FluidIntelligence_0/ukb_step2_FI_0.regenie", header=TRUE)
DB_<-DB[match(c("rs2297086", "rs28496034", "rs75614054"),DB$ID),]
DB_<-DB_%>%mutate(CF="Fluid intelligence")
DBs <- rbind(DBs, DB_)
#Symbol digit substitution
DB <- read.table("SymbolDigit_2/ukb_step2_SDS_2.regenie", header=TRUE)
DB_<-DB[match(c("rs6120880"),DB$ID),]
DB_<-DB_%>%mutate(CF="Symbol digit substitution")
DBs <- rbind(DBs, DB_)

#Vascular Dementia  
DB <- read.table("DemVas/ukb_step2_DemVas.regenie", header=TRUE)
DB_<-DB[match(c("rs2571445"),DB$ID),]
DB_<-DB_%>%mutate(CF="Vascular Dementia")
DBs <- rbind(DBs, DB_)

#P-value
DBs$LOG10P <- as.numeric(DBs$LOG10P)
DBs$P_CF <- 10^(DBs$LOG10P*-1)

DBs
DBs=subset(DBs, select= -c(EXTRA, TEST, CHISQ, LOG10P))
DBs

#Laura's new file containing alleles
DB <- read.table("/data/home/hmy431/outputs/fromLaura/LF_sumstats.txt", header=TRUE)
head(DB)
DB <- subset(DB, select=-c(CHROM, GENPOS, A1FREQ))
DB <- rename(DB, "ALLELE0_LF"="ALLELE0")
DB <- rename(DB, "ALLELE1_LF"="ALLELE1")
DB <- rename(DB, "INFO_LF"="INFO")
str(DB)

#Merging
DBm <- left_join(DBs,DB,by="ID")
dim(DBm)                         #25 21
DBm

#MAF
DBmaf <- read.table("/data/home/hmy431/data/Variants.txt", header=TRUE)
DBmaf <- subset(DBmaf, select=c(ID,MAF))
DBm <- left_join(DBm,DBmaf,by="ID")
DBm

#====Save
saveRDS(DBm, file='/data/home/hmy431/data/Coloclized_Variants.rds')
write.table(DBm, file = "/data/home/hmy431/data/Coloclized_Variants.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write_xlsx(DBm,"/data/home/hmy431/data/Coloclized_Variants.xlsx", col_names = TRUE, format_headers = TRUE)


#Alzheimer's Disease (GWAS)
DB <- read.table("GWAS_AD/ukb_GWASAD_add.txt", sep="\t", header=TRUE)
DB_<-DB[match(c("rs1978487", "rs11865499"),DB$ID),]
DB_<-DB_%>%mutate(CF="Alzheimer's Disease")
DB_<-DB_[-c(2:5)]
DB_<-DB_%>%mutate(INFO=" ",.after=ID)
DBs <- rbind(DBs, DB_)





#============================================Providing complete genetic file including position
setwd("/data/home/hmy431/")
#====Reading Regenie output
DB1 <- read.table("outputs/Dementia/ukb_step2_Dem.regenie", header=TRUE)
dim(DB1)                                #15,640 14
sum(duplicated(DB1$ID))                 #0
#====Reading MAF (UKB-mfi)
DB2 <- read.table("/data/Wolfson-PNU-dementia/lungdev/ukb_maf_v3.txt", header=TRUE)
dim(DB2)                                #93,095,623
sum(duplicated(DB2$ID))                 #320,321
DB2 <- rename(DB2, ALLELE0=Allele2)
DB2 <- rename(DB2, ALLELE1=Allele1)
#====Reading variants with gene name
DB3 <- read.table("outputs/fromLaura/variants_filtered_ID_gene.txt", header=TRUE)
head(DB3)
#====Merging
DB <- merge(DB1, DB2, by=c("ID","ALLELE0","ALLELE1"))
head(DB)
DB <- merge(DB, DB3, by="ID")
head(DB)
#====Prepare
DB=subset(DB, select= c(ID, ALLELE0, ALLELE1, CHROM, GENPOS, MAF, GENE))
DB <- rename(DB, POS=GENPOS)
#====Save
saveRDS(DB, file='data/Variants.rds')
#Tab delimited text file
write.table(DB, file = "data/Variants.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#============================================Investigating mismatch in merging meta-GWAS with UKB
#====Reading downloaded meta-GWAS summary data
library(data.table)
setwd("/data/home/hmy431/")     
df <- fread('/data/home/hmy431/data/GCST90027158_buildGRCh38.tsv.gz')
str(df)                                           #21,101,114 obs. of  17 variables
df <- rename(df, ID=variant_id)
df <- rename(df, position38=base_pair_location)
df <- rename(df, aID=variant_alternate_id)
sum(duplicated(df$position38))                      #1,180,766 position duplicates in GWAS dataset
sum(duplicated(df$ID))                            #30,648 rs No. duplicates in GWAS dataset
names(df)
DBs=subset(df, select= c(ID, chromosome, position38))

#====Reading UKB MAF
DB_MAF <- read.table("/data/Wolfson-PNU-dementia/lungdev/ukb_maf_v3.txt", header=TRUE)
str(DB_MAF)                                       #93,095,623 obs. of  4 variables
sum(duplicated(DB_MAF$ID))      #320,321 Duplicate in UKB-mfi dataset
DB <- merge(df, DB_MAF, by=c("ID","Allele1","Allele2"))
str(DB)                                           #19,067,815

#====Reading variants dataset
DB_V <- read.table("outputs/fromLaura/variants_filtered_ID_gene.txt", header=TRUE)
str(DB_V)                                         #15,640

#====Merging
DBm <- left_join(DB_V,DBs,by="ID")
dim(DBm)                         #15645
sum(is.na(DBm$position38))       #2743

#====Investigating duplicates
sum(duplicated(DBm$ID))                            #5 Duplicate in merged dataset
DBm <- DBm %>% mutate(dup=ifelse(duplicated(DBm$ID)==FALSE, 0,1))
Dup =subset(DBm, dup==1)

#====Clearing duplicates
DBm <- DBm %>% distinct(ID, .keep_all= TRUE)
str(DB)
sum(duplicated(DBm$ID))

#====List of 2743 UKB variants NOT in meta-GWAS
DBmis <- DBm %>% filter(is.na(DBm$position38))
DBmis <- DBmis %>% subset(select=-c(chromosome,position38))
head(DBmis, n=100)

#====Investigating absence of UKB variants in meta-GWAS
# rs556883943   https://www.ncbi.nlm.nih.gov/snp/?term=rs556883943
see<-DB_MAF[match("rs556883943",DB_MAF$ID),]         #In UKB:               YES
see<-DBs[match("rs556883943",DBs$ID),]               #In meta-GWAS, rs:     NA
see<-DBs[match("112524017",DBs$position38),]         #In meta-GWAS, Pos38:  NA    (GRCh38)
see<-DBs[match("113066639",DBs$position38),]         #In meta-GWAS, Pos37:  NA    (GRCh37)
see


#====Save
saveRDS(DBmis, file='data/Missing_ukb_GWASAD.rds')
#Tab delimited text file
write.table(DBmis, file = "data/Missing_ukb_GWASAD.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
