#------------------------------------------------------------------------------
##Colocalization SuSiE:       Code Lab
#------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE))
rm(coloc_FVC, coloc_ratio, DB, DList1, DList2, gDB, Su1, Su2,susie.res)
#Note: Fluid Intelligence as an example
library(dplyr)
library(tidyverse)
library(coloc)
library(writexl)
library(susieR)
library(Tmisc)


#Set working directory
setwd("/data/home/hmy431/")


#LD matrix (once, adding row/column names)      bfile_imp_filt --> zArchive
LDx <- read.table("/data/Wolfson-PNU-dementia/lungdev/zArchive/ld_matrix.ld", header=FALSE)
SNPlist <- read.table("/data/Wolfson-PNU-dementia/lungdev/zArchive/ld_matrix.snplist", header=FALSE)
row.names(LDx) <- SNPlist$V1
colnames(LDx) <- SNPlist$V1
#install.packages("Tmisc")
library(Tmisc)
corner(LDx, n = 5)
LDx <- data.matrix(LDx)
corner(LDx, n = 5)
is.matrix(LDx)
rm(SNPlist)

#LD matrix (once, adding row/column names)      zArchive --> lungdev
LDx1 <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.ld", header=FALSE)
SNPlist1 <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.snplist", header=FALSE)
row.names(LDx1) <- SNPlist1$V1
colnames(LDx1) <- SNPlist1$V1
LDx1 <- data.matrix(LDx1)
corner(LDx1, n = 5)
is.matrix(LDx1)


#write.table(LDx, file = "/data/Wolfson-PNU-dementia/lungdev/ld_matrix.ld", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
#Read LD matrix
#LDx <- read.table("/data/Wolfson-PNU-dementia/lungdev/ld_matrix.ld", header=TRUE)

#Read the input file 
DB <- read.table("outputs/FluidIntelligence_0/ukb_step2_FI_0_add.txt", header=TRUE)
DB <- read.table("outputs/PairsMatching_0/ukb_step2_Pairs_0_add.txt", header=TRUE)
DB <- read.table("outputs/ReactionTime_0/ukb_step2_RT_0_add.txt", header=TRUE)
DB <- read.table("/data/home/hmy431/data/ukb_GWASAD_add.txt", header=TRUE)


#Add position from the varinat file 
DB_P <- read.table("data/Variants.txt", header=TRUE)
head(DB_P)
DB_P <- DB_P[c(1,5)]
DB <- left_join(DB,DB_P,by="ID")
dim(DB)                         #15298
sum(is.na(DB$POS))              #0
rm(DB_P)

#Add position from the bim file 
#DB_P <- read.table("/data/Wolfson-PNU-dementia/lungdev/zArchive/ukb_imp_allChrs_filt.bim", header=FALSE)
#DB_P <- DB_P[c(2,4)]
#DB_P <- rename(DB_P, ID=V2)
#DB_P <- rename(DB_P, POS=V4)
#DB <- left_join(DB,DB_P,by="ID")
#dim(DB)                         #15298
#sum(is.na(DB$POS))              #43
#rm(DB_P)

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
#==colocalisation           Priors: Default   MAF = NULL, p1=1, p2=1e-04, p12=1e-05
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
            
class(coloc_FVC)
coloc_FVC
sensitivity(coloc_FVC,"H4 > 0.9")

tiff("outputs/Figures/FI_FVC_PTCH1.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_FVC,"H4 > 0.9")
dev.off()

#with priors defined         
coloc_FVC2 <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            MAF = NULL, p1=1e-03, p2=1e-04, p12=1e-05) 
class(coloc_FVC2)
coloc_FVC2
sensitivity(coloc_FVC2,"H4 > 0.9")

#with sdY defined
coloc_FVC3 <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, sdY=0.981, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, sdY=0.980 , snp=gDB$ID, position=gDB$POS)) 

class(coloc_FVC3)
coloc_FVC3
sensitivity(coloc_FVC3,"H4 > 0.9")

#==SuSiE                
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=138885, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDx)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Colocalise every pair                    Priors: Default
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2)
  print(susie.res$summary)
}
susie.res$summary

#sensitivity        
tiff("outputs/Figures/FI_FVC_PTCH1_Su.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}
dev.off()

#Colocalise every pair                    Priors: Modified 
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.005667, p2=0.000351, p12=0.000351)
  print(susie.res$summary)
}
susie.res$summary

#sensitivity for SuSiE returns error with modified priors


#==Export to Excel  SuSiE
susum<-susie.res$summary
write_xlsx(list("PTCH1"=susum), "outputs/FluidIntelligence_0/FVC_FI_PTCH1.xlsx", col_names = TRUE, format_headers = TRUE)
#Adding to an existing file
library(openxlsx)
FL <- createWorkbook()
FL <- loadWorkbook("outputs/FluidIntelligence_0/FVC_FI.xlsx")
addWorksheet(FL, genes)
writeData(FL, genes, susum)

saveWorkbook(FL,"outputs/FluidIntelligence_0/FVC_FI.xlsx",overwrite = TRUE)


#==Export to Excel  SuSiE & SNPs
#SuSiE
susum<-susie.res$summary
FL <- createWorkbook()
addWorksheet(FL, "FVC_PTCH1")
writeData(FL, "FVC_PTCH1", susum)
#Summary stat
see <-DB[match(c("rs2227603","rs2227592"),DB$ID),]
see <- see[-c(15:17)]
addWorksheet(FL, "SNPs-FVC-PTCH1")
writeData(FL, "SNPs-FVC-PTCH1", see)

saveWorkbook(FL,"outputs/ReactionTime_0/SuSiE_RT_PTCH1.xlsx",overwrite = TRUE)


#---------------PTCH1: Ratio/FI
#==colocalisation           MAF = NULL, p1=1e-04, p2=1e-04, p12=1e-05
coloc_ratio <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_ratio)
coloc_ratio
tiff("outputs/Figures/RT0_Ratio_NFATC3.tiff", width = 7, height = 5, units = 'in', res = 300, compression = 'lzw')
sensitivity(coloc_ratio,"H4 > 0.9")
dev.off()

#with priors defined    
coloc_ratio2 <-
  coloc.abf(dataset1=list(beta=gDB$beta_ratio, N=306476, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_ratio, varbeta=gDB$se_ratio2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=426102, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS),
            MAF = NULL, p1=1, p2=1e-04, p12=1e-05) 
class(coloc_ratio2)
coloc_ratio2

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

#sensitivity        #Needs positions
if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=DList1,dataset2=DList2)
}




#see variants
see<-DB[match("rs28496034",DB$ID),]
see<-DB[match(c("rs357565", "rs147794976", "rs2297086", "rs75614054","rs28496034"),DB$ID),]
see

write_xlsx(as.data.frame(LDxs), path = "/data/home/hmy431/data/LDx_KAT8_UKB.xlsx", col_names=TRUE)



#=======================LD Matrix
#Extract LD matrix for KAT8 variants
SNPs <- gDB$ID
LDxs <- LDx[SNPs,SNPs]
corner(LDxs, n = 10)
LDxs <- data.matrix(LDxs)
is.matrix(LDxs)
dim(LDxs)

#Extract LD matrix for some variants
SNPs <- c("16:68154994_CA_C", "rs237833", "rs548092276", "16:68154994_CA_C")
LDxs <- LDx[SNPs,SNPs]
corner(LDxs, n = 4)
LDxs <- data.matrix(LDxs)
is.matrix(LDxs)
dim(LDxs)


#Excluding some variants from matrix
gDB <- gDB[-c(1:9),] 
head(gDB)

#Checking LD matrix mismatch (comparing SNP list)
gDB <- gDB  %>% mutate(LDx = ifelse(ID %in% SNPlist$V1, "1", "0"))
table(gDB$LDx)
see<-gDB[match(0,gDB$LDx),]
see                             #rs7188350
see<-SNPlist[match("rs7188350",SNPlist$ID),]
gDB=subset(gDB, LDx==1)                 


#DB <- read.table("outputs/SymbolDigit_2/ukb_step2_SDS_2_add.txt", header=TRUE)
DB4 <- read.table("fromLaura/variants_filtered_ID_gene.txt", header=TRUE)
head(DB4)

genes <- c('MMP24')
gDB <- DB4  %>% mutate(exclusion = ifelse(GENE %in% genes, "1", "0"))
table(gDB$exclusion)
gDB<-filter(gDB, exclusion==1)
gDB<- subset(gDB, select=-c(exclusion))

DBx=subset(gDB, select= c(ID))
head(DBx)
LDx <- ld_matrix(DBx, with_alleles = FALSE, pop = "EUR")


#========External     Note: Did not worked in OnDemand system
install.packages(c("devtools", "knitr", "rmarkdown"))
library(devtools)
install_github(c("MRCIEU/TwoSampleMR","MRCIEU/MRInstruments"))

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
system.file(package='MRCIEU/TwoSampleMR')

#LD Matrix (creates an LD matrix of r values [signed, and not squared])
snps <- c(rs357565, rs147794976, rs2297086, rs75614054)
ieugwasr::ld_reflookup(snp)
library("TwoSampleMR")
ld_matrix(snps, with_alleles = TRUE, pop = "EUR")
#=> The following variants are not present in the LD reference panel rs147794976


#========External: EUR uploaded     Note: created on personal computer
LDxs <- read.table("/data/home/hmy431/data/LDx_KAT8_EUR.Rdata")
LDxs <- data.matrix(LDxs)
is.matrix(LDxs)

#gDB=subset(gDB, ID!="16:31140187_AG_A" && ID!="rs113125305" && ID!="rs138259061")
gDB=subset(gDB, ID!="16:31140187_AG_A")
gDB=subset(gDB, ID!="rs113125305")
gDB=subset(gDB, ID!="rs138259061")


#==colocalisation           
coloc_FVC <-
  coloc.abf(dataset1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS),
            dataset2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS)) 
class(coloc_FVC)
coloc_FVC
sensitivity(coloc_FVC,"H4 > 0.9")


#==SuSiE with defined priors                
DList1=list(beta=gDB$beta_fvc, N=306455, type="quant", MAF=gDB$MAF, pvalue1=gDB$p_fvc, varbeta=gDB$se_fvc2, snp=gDB$ID, position=gDB$POS, LD=LDxs)
check_dataset(DList1,req="LD")
Su1=runsusie(DList1)
summary(Su1)

DList2=list(beta=gDB$BETA, N=428609, type="quant", MAF=gDB$MAF, pvalue1=gDB$pval, varbeta=gDB$SE2, snp=gDB$ID, position=gDB$POS, LD=LDxs)
check_dataset(DList2,req="LD")
Su2=runsusie(DList2)
summary(Su2)

#Modified priors for SuSiE
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(Su1,Su2, p1=0.001728, p2=0.000107, p12=0.000107)
  print(susie.res$summary)
}
susie.res$summary