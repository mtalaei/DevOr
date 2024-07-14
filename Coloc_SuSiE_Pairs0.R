#------------------------------------------------------------------------------
##Colocalization SuSiE:       Pairs matching
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
DB <- read.table("outputs/PairsMatching_0/ukb_step2_Pairs_0_add.txt", header=TRUE)

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




#===============KAT8          Mismatch error at both traits for SuSiE
#Gene 
genes <- c('KAT8')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))
#gDB<-filter(gDB, MAF>=0.05)

#---------------KAT8: FVC/Pairs
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/Pairs0_FVC_KAT8.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.036957	0.002288	0.002288
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.036957, p2=0.002288, p12=0.002288)  
class(coloc_FVC)
coloc_FVC

#SNP Look up
FVC_r.1<-subset(coloc_FVC$results,SNP.PP.H4>0.1)
FVC_r.1                 # rs138259061 (SNP.PP.H4: 0.546), rs1978485 (0.145), rs1978487 (0.146)
see<-DB[match(c("rs138259061", "rs1978485", "rs1978487"),DB$ID),]
see

#==SuSiE                
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "FVC_KAT8")
writeData(FL, "FVC_KAT8", susum)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}


#==Export to Excel
susum<-susie.res$summary
FL <- createWorkbook()
addWorksheet(FL, "Ratio-KAT8")
writeData(FL, "Ratio-KAT8", susum)
saveWorkbook(FL,"outputs/PairsMatching_0/SuSiE_Pairs_KAT8.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}

#---------------KAT8: Ratio/Pairs       NOT applicable 

#=====Export selected summary stat to excel
see <-DB[match(c("", "", "", "",""),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("",""),DB$ID),]
see2 <- see2[-c(15:17)]
addWorksheet(FL, "SNPs-Ratio-KAT8")
addWorksheet(FL, "SNPs-FVC-KAT8")
writeData(FL, "SNPs-Ratio-KAT8", see)
writeData(FL, "SNPs-FVC-KAT8", see2)
saveWorkbook(FL,"outputs/PairsMatching_0/SuSiE_Pairs_KAT8.xlsx",overwrite = TRUE)




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
#see<-SNPlist[match("NA",SNPlist$ID),]
#gDB=subset(gDB, LDx==1)
gDB=subset(gDB, select= -c(LDx))


#---------------CSNK2B: FVC/Pairs       NOT applicable because the gene not in the list for FVC (run exceptionally) 
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/Pairs0_FVC_CSNK2B.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.094444	0.005848	0.005848
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.094444, p2=0.005848, p12=0.005848)  
class(coloc_FVC)
coloc_FVC

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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


#---------------CSNK2B: Ratio/Pairs
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/Pairs0_Ratio_CSNK2B.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.094444	0.005848	0.005848
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.094444, p2=0.005848, p12=0.005848) 
class(coloc_ratio)
coloc_ratio

#==SuSiE                
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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

#sensitivity  
tiff("outputs/Figures/Pairs0_Ratio_CSNK2B_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
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

saveWorkbook(FL,"outputs/PairsMatching_0/SuSiE_Pairs_CSNK2B.xlsx",overwrite = TRUE)



#===============ITGAV             Mismatch error at both traits 
#Gene 
genes <- c('ITGAV')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#Consistency with LD matrix               => confirmed, NO need to run this part
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)
# see<-gDB[match(0,gDB$LDx),]
# see                             #NA
# see<-SNPlist[match("NA",SNPlist$ID),]
# gDB=subset(gDB, LDx==1)
gDB=subset(gDB, select= -c(LDx))

#---------------ITGAV: FVC/Pairs       NOT applicable 

#---------------ITGAV: Ratio/Pairs      Mismatch error (for both traits) for SuSiE
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/Pairs0_Ratio_ITGAV.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.005280	0.000327	0.000327
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.005280, p2=0.000327, p12=0.000327)  
class(coloc_ratio)
coloc_ratio

#SNP Look up
Ratio_r.1<-subset(coloc_ratio$results,SNP.PP.H4>0.1)
Ratio_r.1                 # rs2084448 (SNP.PP.H4: 0.188), rs3816386 (0.167), rs11685758 (0.146)
see<-DB[match(c("rs2084448", "rs3816386", "rs11685758"),DB$ID),]
see

#==SuSiE                
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "Ratio-ITGAV")
writeData(FL, "Ratio-ITGAV", susum)
saveWorkbook(FL,"outputs/PairsMatching_0/SuSiE_Pairs_ITGAV.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}


#=====Export selected summary stat to excel
see <-DB[match(c("", "", "", "",""),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("",""),DB$ID),]
see2 <- see2[-c(15:17)]
#FL <- loadWorkbook("outputs/PairsMatching_0/SuSiE_Pairs_ITGAV.xlsx")
addWorksheet(FL, "SNPs-Ratio-ITGAV")
addWorksheet(FL, "SNPs-FVC-ITGAV")
writeData(FL, "SNPs-Ratio-ITGAV", see)
writeData(FL, "SNPs-FVC-ITGAV", see2)
saveWorkbook(FL,"outputs/PairsMatching_0/SuSiE_Pairs_ITGAV.xlsx",overwrite = TRUE)



#===============================================Summary for KAT8 and ITGAV based on coloc
#=====Export selected summary stat to excel
see <-DB[match(c("rs138259061", "rs1978485", "rs1978487"),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("rs2084448", "rs3816386", "rs11685758"),DB$ID),]
see2 <- see2[-c(15:17)]
FL <- createWorkbook()
addWorksheet(FL, "SNPs-FVC-KAT8")
addWorksheet(FL, "SNPs-Ratio-ITGAV")
writeData(FL, "SNPs-FVC-KAT8", see)
writeData(FL, "SNPs-Ratio-ITGAV", see2)
saveWorkbook(FL,"outputs/PairsMatching_0/SNPs_Pairs_KAT8_ITGAV.xlsx",overwrite = TRUE)