#------------------------------------------------------------------------------
##Colocalization SuSiE:       Reaction time 
#------------------------------------------------------------------------------

rm(list = ls(all.names = TRUE))
rm(coloc_FVC, coloc_ratio, DList1, DList2, gDB, Su1, Su2, susie.res, susum, see, see2, FL)

library(dplyr)
library(tidyverse)
library(coloc)
#library(writexl)
library(openxlsx)
#library(susieR)


#Set working directory
setwd("/data/home/hmy431/")


LDx <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.ld", header=FALSE)
SNPlist <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.snplist", header=FALSE)
row.names(LDx) <- SNPlist$V1
colnames(LDx) <- SNPlist$V1
LDx <- data.matrix(LDx)
is.matrix(LDx)
dim(LDx)
rm(SNPlist)


#Read the input file 
DB <- read.table("outputs/ReactionTime_0/ukb_step2_RT_0_add.txt", header=TRUE)

#Add position from the variant file 
DB_P <- read.table("data/Variants.txt", header=TRUE)
head(DB_P)
DB_P <- DB_P[c(1,5)]
DB <- left_join(DB,DB_P,by="ID")
dim(DB)                         #15298
sum(is.na(DB$POS))              #0
rm(DB_P)

#Varbeta
DB$SE2 <- DB$SE^2 
DB$se_fvc2 <- DB$se_fvc^2 
DB$se_ratio2 <- DB$se_ratio^2 
sum(is.na(DB$beta_fvc))       #2 missings (2 rows)
DB <- na.omit(DB)             #no missing 
sum(is.na(DB$beta_ratio))     #0 missings




#===============PTCH1
#Gene 
genes <- c('PTCH1')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))


#---------------PTCH1: FVC/RT
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/RT0_FVC_PTCH1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.005667	0.000351	0.000351
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.005667, p2=0.000351, p12=0.000351) 
class(coloc_FVC)
coloc_FVC

#==SuSiE                Priors: Default
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Colocalise every pair
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2)
  print(susie.res$summary)
}
susie.res$summary

#==Export to Excel
susum<-susie.res$summary
FL <- createWorkbook()
addWorksheet(FL, "FVC_PTCH1")
writeData(FL, "FVC_PTCH1", susum)

#sensitivity        
tiff("outputs/Figures/RT0_FVC_PTCH1_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()

#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.005667, p2=0.000351, p12=0.000351)
  print(susie.res$summary)
}
susie.res$summary


#---------------PTCH1: Ratio/RT
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/RT0_Ratio_PTCH1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.005667	0.000351	0.000351
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.005667, p2=0.000351, p12=0.000351) 
class(coloc_ratio)
coloc_ratio

FVC_r.1<-subset(coloc_ratio$results,SNP.PP.H4>0.1)
FVC_r.1

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Colocalise every pair
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2)
  print(susie.res$summary)
}
susie.res$summary

#==Export to Excel
susum<-susie.res$summary
#FL <- createWorkbook()
addWorksheet(FL, "Ratio-PTCH1")
writeData(FL, "Ratio-PTCH1", susum)
#saveWorkbook(FL,"outputs/ReactionTime_0/SuSiE_RT_PTCH1.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
}
tiff("outputs/Figures/RT0_Ratio_PTCH1_SuR1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/RT0_Ratio_PTCH1_SuR2.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
dev.off()

#=====Export selected summary stat to excel
see <-DB[match(c("rs113154802","rs539078574"),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("rs2297086", "rs539078574","rs75614054"),DB$ID),]
see2 <- see2[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-PTCH1")
addWorksheet(FL, "SNPs-Ratio-PTCH1")
writeData(FL, "SNPs-FVC-PTCH1", see)
writeData(FL, "SNPs-Ratio-PTCH1", see2)
saveWorkbook(FL,"outputs/ReactionTime_0/SuSiE_RT_PTCH1.xlsx",overwrite = TRUE)


#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.005667, p2=0.000351, p12=0.000351)
  print(susie.res$summary)
}
susie.res$summary


#===============SERPINC1
#Gene 
genes <- c('SERPINC1')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))


#---------------SERPINC1: FVC/RT
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
#tiff("outputs/Figures/RT0_FVC_SERPINC1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
#sensitivity(coloc_FVC,"H4 > 0.9")
#dev.off()

# 0.026563	0.001645	0.001645
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.026563, p2=0.001645, p12=0.001645) 
class(coloc_FVC)
coloc_FVC

#==SuSiE
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Colocalise every pair                Priors: Default
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2)
  print(susie.res$summary)
}
susie.res$summary

#sensitivity        
tiff("outputs/Figures/RT0_FVC_SERPINC1_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()


#Colocalise every pair                Priors:Modified 0.026563	0.001645	0.001645
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.026563, p2=0.001645, p12=0.001645)
  print(susie.res$summary)
}
susie.res$summary

#==Export to Excel
#SuSiE
susum<-susie.res$summary
FL <- createWorkbook()
addWorksheet(FL, "FVC_SERPINC1")
writeData(FL, "FVC_SERPINC1", susum)
#Summary stat
see <-DB[match(c("rs2227603","rs2227592"),DB$ID),]
see <- see[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-SERPINC1")
writeData(FL, "SNPs-FVC-SERPINC1", see)

saveWorkbook(FL,"outputs/ReactionTime_0/SuSiE_RT_SERPINC1_m.xlsx",overwrite = TRUE)
