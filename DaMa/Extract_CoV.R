##Extracting co-variables: age, sex, batch, center, and PC variables  
library(tidyverse)
library(readr)
library(dplyr)
library(reshape)

#All variables at instance 2 
df <- read_tsv("/data/Wolfson-PNU-dementia/datasets/ukb_extracted/47532_Coded/1_ukbpheno_CoreDataClin1.tsv") 
covDB = df%>% select(contains("EID"),
                     contains("Age when attended assessment centre.0.0"),
                     contains("Sex.0.0"), 
                     contains("UK Biobank assessment centre.0.0"),
                     contains("Genotype measurement batch.0.0"), 
                     contains("Genetic principal components"))
covDB <- rename(covDB, c(EID="IID"))
covDB$FID <- covDB$IID
covDB <- covDB[-c(4,17:47)]
covDB <- covDB %>% relocate(FID, .before = `Age when attended assessment centre.0.0`)
covDB <- covDB %>% relocate(IID, .after =  FID)
#Report
print("Dataset dimension:")
dim(covDB)
print("Dataset structure:")
str(covDB)
print("====Summary Statistics:")
print("Age at baseline")
summary(covDB$`Age when attended assessment centre.0.0`)
print("Sex")
table(covDB$Sex.0.0)
print("Batch")
table(covDB$`Genotype measurement batch.0.0`)
print("Assessment centre")
table(covDB$`UK Biobank assessment centre.0.0`)
#Rename
covDB <- rename(covDB, c(`Age when attended assessment centre.0.0`="age",
                         `Sex.0.0`="sex",
                         `UK Biobank assessment centre.0.0`="centre",
                         `Genotype measurement batch.0.0`="batch",
                         `Genetic principal components.0.1`="pc1",
                         `Genetic principal components.0.2`="pc2",
                         `Genetic principal components.0.3`="pc3",
                         `Genetic principal components.0.4`="pc4",
                         `Genetic principal components.0.5`="pc5",
                         `Genetic principal components.0.6`="pc6",
                         `Genetic principal components.0.7`="pc7",
                         `Genetic principal components.0.8`="pc8",
                         `Genetic principal components.0.9`="pc9",
                         `Genetic principal components.0.10`="pc10"))
#New variables
covDB <- covDB %>% mutate(array = "") %>% 
mutate(array = ifelse(batch > 22, 3, array)) %>% 
mutate(array = ifelse(batch > 0 & batch <= 22, 2, array)) %>% 
mutate(array = ifelse(batch <= 0, 1, array))
covDB <- covDB %>% relocate(array, .after = sex)
print("Array")
table(covDB$array)
covDB=subset(covDB, select= -c(batch))
#Report
print("Renamed dataset structure:")
str(covDB)
print("Variable names:")
names(covDB)
#Save
saveRDS(covDB, file="/data/Wolfson-PNU-dementia/lungdev/ukb_cov.rds")
print("Saved as: ukb_cov.rds")
write.table(covDB, file = "/data/Wolfson-PNU-dementia/lungdev/ukb_cov.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Saved as: ukb_cov.txt")
