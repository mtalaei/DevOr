#------------------------------------------------------------------------------
##Colocalization SuSiE:       Fluid Intelligence
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


#LD matrix 
LDx <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.ld", header=FALSE)
SNPlist <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.snplist", header=FALSE)
row.names(LDx) <- SNPlist$V1
colnames(LDx) <- SNPlist$V1
LDx <- data.matrix(LDx)
is.matrix(LDx)
dim(LDx)
#rm(SNPlist)


#Read the input file 
DB <- read.table("outputs/FluidIntelligence_0/ukb_step2_FI_0_add.txt", header=TRUE)

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


#---------------PTCH1: FVC/FI
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/FI0_FVC_PTCH1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.005667	0.000351	0.000351
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.005667, p2=0.000351, p12=0.000351) 
class(coloc_FVC)
coloc_FVC

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
tiff("outputs/Figures/FI0_FVC_PTCH1_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
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


#---------------PTCH1: Ratio/FI
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/FI0_Ratio_PTCH1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.005667	0.000351	0.000351
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.005667, p2=0.000351, p12=0.000351) 
class(coloc_ratio)
coloc_ratio

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "Ratio-PTCH1")
writeData(FL, "Ratio-PTCH1", susum)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_PTCH1.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
}
tiff("outputs/Figures/FI0_Ratio_PTCH1_SuR1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_PTCH1_SuR2.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
dev.off()


#=====Export selected summary stat to excel
see <-DB[match(c("rs113154802","rs28496034"),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("rs2297086", "rs75614054","rs28496034"),DB$ID),]
see2 <- see2[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-PTCH1")
addWorksheet(FL, "SNPs-Ratio-PTCH1")
writeData(FL, "SNPs-FVC-PTCH1", see)
writeData(FL, "SNPs-Ratio-PTCH1", see2)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_PTCH1.xlsx",overwrite = TRUE)


#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.005667, p2=0.000351, p12=0.000351)
  print(susie.res$summary)
}
susie.res$summary


#===============CSNK2B 
rm(coloc_FVC, coloc_ratio, DList1, DList2, gDB, Su1, Su2, susie.res, susum, see, see2, FL)

#Gene 
genes <- c('CSNK2B')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#Consistency with LD matrix               => confirmed, NO need to run this part
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)                  #1=>9
see<-gDB[match(0,gDB$LDx),]
see                             #NA
see<-SNPlist[match("NA",SNPlist$ID),]
gDB=subset(gDB, LDx==1)
gDB=subset(gDB, select= -c(LDx))


#---------------CSNK2B: FVC/FI       NOT applicable because the gene not in the list for FVC (run exceptionally)
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/FI0_FVC_CSNK2B.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.094444	0.005848	0.005848
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.094444, p2=0.005848, p12=0.005848) 
class(coloc_FVC)
coloc_FVC

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "FVC_CSNK2B")
writeData(FL, "FVC_CSNK2B", susum)

#sensitivity  
tiff("outputs/Figures/FI0_FVC_CSNK2B_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()

#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.094444, p2=0.005848, p12=0.005848)
  print(susie.res$summary)
}
susie.res$summary

see<-DB[match(c("rs9267531"),DB$ID),]
see 

#==Export to Excel
susum<-susie.res$summary
addWorksheet(FL, "FVC_CSNK2B_m")
writeData(FL, "FVC_CSNK2B_m", susum)


#---------------CSNK2B: Ratio/FI
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/FI0_Ratio_CSNK2B.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.094444	0.005848	0.005848
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.094444, p2=0.005848, p12=0.005848) 
class(coloc_ratio)
coloc_ratio

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "Ratio-CSNK2B")
writeData(FL, "Ratio-CSNK2B", susum)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_CSNK2B.xlsx",overwrite = TRUE)

#sensitivity  
tiff("outputs/Figures/FI0_Ratio_CSNK2B_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()


#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.094444, p2=0.005848, p12=0.005848)
  print(susie.res$summary)
}
susie.res$summary

#==Export to Excel
susum<-susie.res$summary
addWorksheet(FL, "Ratio-CSNK2B_m")
writeData(FL, "Ratio-CSNK2B_m", susum)

#=====Export selected summary stat to excel
see <-DB[match(c("rs9267531","rs3117578","rs575597758"),DB$ID),]
see <- see[-c(15:17)]
addWorksheet(FL, "SNPs-CSNK2B")
writeData(FL, "SNPs-CSNK2B", see)

saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_CSNK2B.xlsx",overwrite = TRUE)



#===============NFATC3 UNDER Suspicion 
rm(coloc_FVC, coloc_ratio, DList1, DList2, gDB, Su1, Su2, susie.res, susum, see, see2, FL)

#Gene 
genes <- c('NFATC3')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#Consistency with LD matrix
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)
see<-gDB[match(0,gDB$LDx),]
see                             #rs7188350
see<-SNPlist[match("rs7188350",SNPlist$ID),]
gDB=subset(gDB, LDx==1)
gDB=subset(gDB, select= -c(LDx))


#---------------NFATC3: FVC/FI       NOT applicable because the gene not in the list for FVC (run exceptionally)
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/FI0_FVC_NFATC3.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.003025	0.000187	0.000187
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.003025, p2=0.000187, p12=0.000187) 
class(coloc_FVC)
coloc_FVC

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "FVC_NFATC3")
writeData(FL, "FVC_NFATC3", susum)

#sensitivity  
tiff("outputs/Figures/FI0_FVC_NFATC3_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()

#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.003025, p2=0.000187, p12=0.000187)
  print(susie.res$summary)
}
susie.res$summary

see<-DB[match(c("rs237833"),DB$ID),]
see 


#---------------NFATC3: Ratio/FI
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/FI0_Ratio_NFATC3.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.003025	0.000187	0.000187
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.003025, p2=0.000187, p12=0.000187) 
class(coloc_ratio)
coloc_ratio

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "Ratio-NFATC3")
writeData(FL, "Ratio-NFATC3", susum)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_NFATC3.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=3,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=4,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=5,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=6,dataset1=DList1,dataset2=DList2)
}
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR2.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR3.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=3,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR4.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=4,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR5.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=5,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/FI0_Ratio_NFATC3_SuR6.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=6,dataset1=DList1,dataset2=DList2)
dev.off()

#=====Export selected summary stat to excel
see <-gDB[match(c("16:68154994_CA_C", "rs237833", "rs548092276", "16:68154994_CA_C", "rs11275011"),gDB$ID),]
see <- see[-c(15:17)]
addWorksheet(FL, "SNPs-Ratio-NFATC3")
writeData(FL, "SNPs-Ratio-NFATC3", see)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_NFATC3.xlsx",overwrite = TRUE)


#Modified priors  for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.003025, p2=0.000187, p12=0.000187)
  print(susie.res$summary)
}
susie.res$summary



#Extract LD matrix (Because only rs237833 is recognized-rs548092276 rs11275011 is missing from 1000G (GRCh37))
SNPs <- c("16:68154994_CA_C", "rs237833", "rs548092276", "16:68154994_CA_C", "rs11275011")
LDxs <- LDx[SNPs,SNPs]
corner(LDxs, n = 5)
LDxs <- data.matrix(LDxs)
is.matrix(LDxs)
dim(LDxs)