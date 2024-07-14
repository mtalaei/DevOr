#------------------------------------------------------------------------------
##Colocalization SuSiE:       Meta-GWAS Alzheimer's Disease 
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


#LD matrix (once, adding row/column names)      bfile_imp_filt --> zArchive
LDx <- read.table("/data/Wolfson-PNU-dementia/lungdev/zArchive/ld_matrix.ld", header=FALSE)
SNPlist <- read.table("/data/Wolfson-PNU-dementia/lungdev/zArchive/ld_matrix.snplist", header=FALSE)
row.names(LDx) <- SNPlist$V1
colnames(LDx) <- SNPlist$V1
LDx <- data.matrix(LDx)
is.matrix(LDx)
dim(LDx)
rm(SNPlist)


#Read the input file 
DB <- read.table("/data/home/hmy431/data/ukb_GWASAD_add.txt", header=TRUE)

#Add position from the variant file 
DB_P <- read.table("data/Variants.txt", header=TRUE)
head(DB_P)
DB_P <- DB_P[c(1,5)]
DB <- left_join(DB,DB_P,by="ID")
dim(DB)                         #15297
sum(is.na(DB$POS))              #0
rm(DB_P)

#Varbeta
DB$SE2 <- DB$SE^2 
DB$se_fvc2 <- DB$se_fvc^2 
DB$se_ratio2 <- DB$se_ratio^2 
sum(is.na(DB$beta_fvc))       #2 missings (2 rows)
DB <- na.omit(DB)             #no missing 
sum(is.na(DB$beta_ratio))     #0 missings


#===============KAT8
#Gene 
genes <- c('KAT8')
gDB <- DB  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

#Consistency with LD matrix
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)                  #1=>16
see<-gDB[match(0,gDB$LDx),]
see                             #NA
gDB=subset(gDB, select= -c(LDx))


#---------------KAT8: FVC/AD            Mismatch error (for both traits) for SuSiE
#==colocalisation           MAF = NULL, p1=1, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=788989, s=0.14109, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
tiff("outputs/Figures/AD_FVC_KAT8.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

# P1=P2=1e-04 p12=5e-06 
coloc_FVC2 <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID),
            dataset2=list(beta=gDB$BETA, N=788989, s=0.14109, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID),
            MAF = NULL, p1=1e-04, p2=1e-04, p12=5e-06) 
class(coloc_FVC2)
coloc_FVC2
sensitivity(coloc_FVC2,"H4 > 0.9")

# 0.053125   0.0032895   0.0032895 (KAT8 with 16 SNPs rather than original 23) 
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID),
            dataset2=list(beta=gDB$BETA, N=788989, s=0.14109, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID),
            MAF = NULL, p1=0.053125, p2=0.0032895, p12=0.0032895) 
class(coloc_FVC)
coloc_FVC

#SNP Look up
FVC_r.1<-subset(coloc_FVC$results,SNP.PP.H4>0.1)
FVC_r.1                 # rs1978487 (SNP.PP.H4: 0.318), rs1978485 (0.261; R2 0.97), rs1549293 (0.104; R2 0.93)
see<-DB[match(c("rs1978487", "rs1978485", "rs1549293"),DB$ID),]

#==SuSiE                without priors defined
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=788989, s=0.14109, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
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

#sensitivity        #Needs positions
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}


#==Fine mapping: FVC    with modified priors p1=0.053125, p2=0.0032895
#(KAT8 with 16 SNPs rather than original 23)
fmap_FVC <-
  finemap.abf(dataset=list(beta=gDB$beta_fvc, varbeta=gDB$se_fvc2, N=306455, 
                           type="quant", MAF=gDB$MAF, pvalues=gDB$p_fvc, snp=gDB$ID), p1=0.053125)
fmap_FVC <- na.omit(fmap_FVC)
fmap_FVC=subset(fmap_FVC, select= c(snp,SNP.PP))
fmap_FVC <- rename(fmap_FVC, SNP.PP.FVC=SNP.PP)
str(fmap_FVC)
#preparing association data
gDBs=subset(gDB, select= c(ID,p_fvc,beta_fvc,se_fvc2))
gDBs <- rename(gDBs, snp=ID)
#merge with association data
fmap_FVC <- merge(gDBs,fmap_FVC,by="snp")
str(fmap_FVC)

#==Fine mapping: AD
fmap_AD <-
  finemap.abf(dataset=list(beta=gDB$BETA, varbeta=gDB$SE2, N=487511, s=0.2139, 
                           type="cc", MAF=gDB$MAF, pvalues=gDB$pval, snp=gDB$ID), p1=0.0032895)
fmap_AD <- na.omit(fmap_AD)
fmap_AD=subset(fmap_AD, select= c(snp,SNP.PP))
fmap_AD <- rename(fmap_AD, SNP.PP.AD=SNP.PP)
#preparing association data
gDBs=subset(gDB, select= c(ID,pval,BETA,SE2))
gDBs <- rename(gDBs, snp=ID)
#merge with association data
fmap_AD <- merge(gDBs,fmap_AD,by="snp")
str(fmap_AD)

#----merging FVC & AD
fmap_FVC_AD <- merge(fmap_FVC,fmap_AD,by="snp")
str(fmap_FVC_AD)

#selection PP>10%
fmap_FVC_AD.1 <- subset(fmap_FVC_AD,SNP.PP.FVC>0.1 | SNP.PP.AD>0.1)
fmap_FVC_AD.1


#---------------KAT8: Ratio/FI       NOT applicable (KAT8 is not linked to the ratio)
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=788989, s=0.14109, type="cc", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/AD_Ratio_KAT8.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()


#==Export selected summary stat to excel          NOT UPDATED FOR AD-KAT8
see <-DB[match(c("", "", "", "",""),DB$ID),]
see <- see[-c(15:17)]
see2 <-DB[match(c("",""),DB$ID),]
see2 <- see2[-c(15:17)]
#FL <- loadWorkbook("outputs/FluidIntelligence_0/SuSiE_FI_KAT8.xlsx")
addWorksheet(FL, "SNPs-Ratio-KAT8")
addWorksheet(FL, "SNPs-FVC-KAT8")
writeData(FL, "SNPs-Ratio-KAT8", see)
writeData(FL, "SNPs-FVC-KAT8", see2)
saveWorkbook(FL,"outputs/FluidIntelligence_0/SuSiE_FI_KAT8.xlsx",overwrite = TRUE)




