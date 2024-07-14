#------------------------------------------------------------------------------
##Colocalization SuSiE:       Vascular Dementia
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
#rm(SNPlist)


#Read the input file 
DB <- read.table("outputs/DemVas/ukb_step2_DemVas_add.txt", header=TRUE)


#Add position from the variant file 
#DB_P <- read.table("data/Variants.txt", header=TRUE)
#head(DB_P)
#DB_P <- DB_P[c(1,5)]
#DB <- left_join(DB,DB_P,by="ID")
#dim(DB)                         #15298
#sum(is.na(DB$POS))              #0
#rm(DB_P)

#Varbeta
DB$SE2 <- DB$SE^2 
DB$se_fvc2 <- DB$se_fvc^2 
DB$se_ratio2 <- DB$se_ratio^2 
sum(is.na(DB$beta_fvc))       #2 missings (2 rows)
DB <- na.omit(DB)             #no missing 
sum(is.na(DB$beta_ratio))     #0 missings




#===============TNS1
#Gene 
genes <- c('TNS1')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#Consistency with LD matrix
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)                  #1/491
see<-gDB[match(0,gDB$LDx),]
see                             #rs34572447
#see<-SNPlist[match("NA",SNPlist$ID),]
gDB=subset(gDB, LDx==1)
gDB=subset(gDB, select= -c(LDx))


#---------------TNS1: FVC/DemVas
#Note: Comprable values for H1 and H4: 49% and 48%, respectively 
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/DemVas_FVC_TNS1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.001728	0.000107	0.000107
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.001728, p2=0.000107, p12=0.000107) 
class(coloc_FVC)
coloc_FVC

#==SuSiE                Priors: Default               NOT applicable; NOT DONE till Ratio
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "FVC_TNS1")
writeData(FL, "FVC_TNS1", susum)

#sensitivity        
tiff("outputs/Figures/DemVas_FVC_TNS1_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()


#---------------TNS1: Ratio/DemVas
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/DemVas_Ratio_TNS1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

# 0.001728	0.000107	0.000107
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.001728, p2=0.000107, p12=0.000107) 
class(coloc_ratio)
coloc_ratio

#SNP Look up
Ratio_r.1<-subset(coloc_ratio$results,SNP.PP.H4>0.1)
Ratio_r.1                 # rs2571445 (SNP.PP.H4: 0.470)
see<-DB[match(c("rs2571445"),DB$ID),]
see

#==SuSiE                Mismatch error for the Ratio but NO credible set for VasDem          
DList1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=422785, s=0.00365, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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
addWorksheet(FL, "Ratio-TNS1")
writeData(FL, "Ratio-TNS1", susum)
#saveWorkbook(FL,"outputs/DemVas/SuSiE_RT_TNS1.xlsx",overwrite = TRUE)

#sensitivity        
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
}
tiff("outputs/Figures/DemVas_Ratio_TNS1_SuR1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
dev.off()
tiff("outputs/Figures/DemVas_Ratio_TNS1_SuR2.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=DList1,dataset2=DList2)
dev.off()

#=====Export selected summary stat to excel
see <-DB[match(c("",""),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("", "",""),DB$ID),]
see2 <- see2[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-TNS1")
addWorksheet(FL, "SNPs-Ratio-TNS1")
writeData(FL, "SNPs-FVC-TNS1", see)
writeData(FL, "SNPs-Ratio-TNS1", see2)
saveWorkbook(FL,"outputs/DemVas/SuSiE_RT_TNS1.xlsx",overwrite = TRUE)


