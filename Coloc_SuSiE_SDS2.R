#------------------------------------------------------------------------------
##Colocalization SuSiE:       Symbol digit substitution
#------------------------------------------------------------------------------

rm(list = ls(all.names = TRUE))
rm(coloc_FVC, coloc_ratio, DB, DList1, DList2, gDB, Su1, Su2,susie.res)

library(dplyr)
library(tidyverse)
library(coloc)
#library(writexl)
library(openxlsx)
#library(susieR)


#Set working directory
setwd("/data/home/hmy431/")


#LD matrix for MMP24 (Provided by Laura)      
LDx <- read.table("/data/home/hmy431/data/mmp24_ld.txt", header=TRUE)
library(Tmisc)
corner(LDx, n = 8)
colnames(LDx) <- row.names(LDx)
corner(LDx, n = 8)
LDx <- data.matrix(LDx)
is.matrix(LDx)
dim(LDx)


#Read the input file 
DB <- read.table("outputs/SymbolDigit_2/ukb_step2_SDS_2_add.txt", header=TRUE)

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




#===============MMP24       
#Gene 
genes <- c('MMP24')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#---------------MMP24: FVC/SDS
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=33192, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/SDS2_FVC_MMP24.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# 0.009444	0.000585	0.000585
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=33192, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            p1=0.009444, p2=0.000585, p12=0.000585) 
class(coloc_FVC)
coloc_FVC

#SNP Look up
FVC_r.1<-subset(coloc_FVC$results,SNP.PP.H4>0.1)
FVC_r.1                 # rs6120880 (SNP.PP.H4: 0.986)
see<-DB[match(c("rs6120880"),DB$ID),]
see

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=33192, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Colocalise every pair
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2)
  print(susie.res$summary)
}
susie.res$summary

#==Export to Excel                    NOT Applicable
susum<-susie.res$summary
FL <- createWorkbook()
addWorksheet(FL, "FVC_MMP24")
writeData(FL, "FVC_MMP24", susum)

#sensitivity                          NOT Applicable
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}


#---------------MMP24: Ratio/SDS       NOT Applicable
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=33192, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/SDS2_Ratio_MMP24.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

#==SuSiE                



#=====Export selected summary stat to excel       NOT Applicable
see <-DB[match(c("rs6120880"),DB$ID),]
see <- see[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-MMP24")
writeData(FL, "SNPs-FVC-MMP24", see)
saveWorkbook(FL,"outputs/SymbolDigit_2/SuSiE_SDS_MMP24.xlsx",overwrite = TRUE)

