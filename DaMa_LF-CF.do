*===========================================================================
* Project: Shared Developmental Pathway (UKB)
* Purpose: Models for LuFu-CogFu associations
* Output files: Tables_LF-CF.xlsx
* Launch date: 08 Apr 2022 
* Last update: 04 June 2024
*===========================================================================

/*QMUL*/

/*Office*/ cd "C:\Users\hmy431\OneDrive - Queen Mary, University of London\DevOr\"
/*Home*/ cd "C:\Users\mtala\OneDrive - Queen Mary, University of London\DevOr\"
use "data\ukb_LF0-DmCF.dta", clear
use "data\ukb_LF0-DmCF_m.dta", clear


*===========Compress
set more off
compress


*=================================Exposures
count if FEV1FVC_b0>=1 & !missing(FEV1FVC_b0)			/*n=1*/
*===========FEV1 
codebook FEV1_b0
replace FEV1_b0=. if FEV1FVC_b0>=1
sum FEV1_b0, d 
la var FEV1_b0 "FEV1 (L) best measure @Baseline"
hist FEV1_b0

egen FEV1_b0q5=cut(FEV1_b0), group(5) icode
la var FEV1_b0q5 "Quintiles of FEV1 best measure @Baseline"
tab FEV1_b0q5
tabstat FEV1_b0, by(FEV1_b0q5) statistic(min max) long format

tabstat FEV1_b0, by(FEV1_b0q5) statistic(median mean n min max) format(%9.2f)
recode FEV1_b0q5 (0=1.93) (1=2.39) (2=2.77) (3=3.21) (4=3.91), gen(FEV1_b0q5m)
tab FEV1_b0q5 FEV1_b0q5m

*===========FVC 
codebook FVC_b0
replace FVC_b0=. if FEV1FVC_b0>=1
sum FVC_b0, d 
la var FVC_b0 "FVC (L) best measure @Baseline"
hist FVC_b0

egen FVC_b0q5=cut(FVC_b0), group(5) icode
la var FVC_b0q5 "Quintiles of FVC best measure @Baseline"
tab FVC_b0q5
tabstat FVC_b0, by(FVC_b0q5) statistic(min max) long format

tabstat FVC_b0, by(FVC_b0q5) statistic(median mean n min max) format(%9.2f)
recode FVC_b0q5 (0=2.62) (1=3.17) (2=3.65) (3=4.24) (4=5.12), gen(FVC_b0q5m)
tab FVC_b0q5 FVC_b0q5m

*===========FEV1/FVC 
codebook FEV1FVC_b0
replace FEV1FVC_b0=. if FEV1FVC_b0>=1
sum FEV1FVC_b0, d 
la var FEV1FVC_b0 "FEV1/FVC ratio best measure @Baseline"
hist FEV1FVC_b0

egen FEV1FVC_b0q5=cut(FEV1FVC_b0), group(5) icode
la var FEV1FVC_b0q5 "Quintiles of FEV1/FVC ratio best measure @Baseline"
tab FEV1FVC_b0q5
tabstat FEV1FVC_b0, by(FEV1FVC_b0q5) statistic(min max) long format

tabstat FEV1FVC_b0, by(FEV1FVC_b0q5) statistic(median mean n min max) format(%9.3f)
recode FEV1FVC_b0q5 (0=0.675) (1=0.733) (2=0.764) (3=0.791) (4=0.824), gen(FEV1FVC_b0q5m)
tab FEV1FVC_b0q5 FEV1FVC_b0q5m


*=================================covariates
sum age
tab sex

*Center
codebook centre
tab centre
gen centre_m = centre
recode centre_m 10003=11024 
tab centre_m

*country
*https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=11002
gen country =2 if centre_m==11003 | centre_m==11022 | centre_m==11023 
recode country .=3 if centre_m==11004 | centre_m==11005
recode country .=1

la def country 1"England" 2"Wales" 3"Scotland"
la val country country
tab country

*Height
sum Height_0
hist Height_0
count if Height_0<130							/*19*/
hist Height_0 if Height_0>=130

*Ethnicity
codebook Ethn
recode Ethn 1=0 1001/1003=0 2=3 2001/2004=3 3=1 3001/3004=1 4=2 4001/4003=2 5=3 6=3 -3/-1=4, gen(Ethn_m)
la def Ethn 0"White" 1"South Asian" 2"Black" 3"Mixed, Chinese, other" 4"Unknown"
la val Ethn_m Ethn
tab Ethn_m 
tab Ethn Ethn_m
count if Ethn_m !=0 & !missing(Dementia)			/*0*/

codebook TDI
la var TDI "Townsend deprivation index at recruitment"
sum TDI

/*Edu
1:College or University degree			2:A levels/AS levels or equivalent
3:O levels/GCSEs or equivalent			4:CSEs or equivalent
5:NVQ or HND or HNC or equivalent		6:Other professional eg: nursing, teaching
-7 :None of the above					-3 :Prefer not to answer	*/
tab Edu_0, missing
recode Edu_0 1=0 2=1 3=2 4=3 5=4 6=5 -7=6 -3=6 .=7, gen(Edu_0m)
la var Edu_0m "Qualifications" 
la def Edu 0"College or University degree" 1"A levels/AS levels" 2"O levels/GCSEs" ///
3"CSEs" 4"NVQ or HND or HNC" 5"Other professional qualifications" 6"None of the above" 7"Missing"
la val Edu_0m Edu
tab Edu_0 Edu_0m
tab Edu_0m

recode Edu_0m 3/5=3 6=4 7=5, gen(Edu_0m5)
la var Edu_0m5 "Qualifications" 
la def Edu5 0"College or University degree" 1"A levels/AS levels" 2"O levels/GCSEs" ///
3"CSEs/NVQ/HND/HNC or Other professional qualifications" 4"None of the above" 5"Missing"
la val Edu_0m5 Edu5
tab Edu_0m Edu_0m5
tab Edu_0m5

*Smoking
codebook Smoking_0
recode Smoking_0 (-3=.) (0=0 Never) (1=1 Previous) (2=2 Current), gen(Smoking_0m) test
la var Smoking_0m "Smoking status @baseline"
tab Smoking_0m, missing

*Pack-year
codebook Packyears_0
sum Packyears_0
tabstat Packyears_0, by(Smoking_0m) statistic(n mean median min max) format(%9.2f)

gen Packyears_0m = Packyears_0
recode Packyears_0m .=0 if Smoking_0==0
sum Packyears_0m
tabstat Packyears_0m, by(Smoking_0m) statistic(n mean median min max) format(%9.2f)

count if missing(Packyears_0) & Smoking_0m==1 					/*64,328*/
count if missing(Packyears_0) & Smoking_0m==2 					/*10,774*/
tabstat Packyears_0m, by(Smoking_0m) statistic(n mean median min max) format(%9.2f)

recode Packyears_0m .=17 if Smoking_0m==1
recode Packyears_0m .=24.7 if Smoking_0m==2
sum Packyears_0m

*HTN
la var HTN_0 "HTN @Baseline"
la def ny 0"No" 1"Yes" 
la val HTN_0 ny
tab HTN_0

tab HTN_0 DemVas, row
tab HTN_0 DemAlz, row
mhodds DemVas HTN_0				/*3.42*/
cs DemVas HTN_0
mhodds DemAlz HTN_0				/*2.25*/
cs DemAlz HTN_0


*=================================Merge adjusted Tenure variable
merge 1:1 IID using "data\DementiaTenureA.dta"					/*75 not matched from master*/
replace DementiaTenureA=. if missing(Dementia)
count if missing(DementiaTenureA) & !missing(DementiaTenure)	/*18*/
replace DementiaTenureA=DementiaTenure if missing(DementiaTenureA) & !missing(DementiaTenure)
sum DementiaTenure DementiaTenureA
la var DementiaTenure "Years from baseline to dementia diagnosis"
la var DementiaTenureA "Years from baseline to dementia diagnosis, modified"

count if DementiaTenure<=0 & Dementia_Incident==1				/*0*/		
count if DementiaTenureA<=0 & Dementia_Incident==1				/*0*/		
*list Dementia pDementia_0 iDementia DementiaTenure DementiaTenureA if DementiaTenureA<=0 & !missing(iDementia)
*replace DementiaTenureA=. if DementiaTenureA<0
drop _merge

*Cleaning
count if DementiaTenureA <0
gen DementiaTenureAA=DementiaTenureA
replace DementiaTenureAA=. if DementiaTenureAA<0
/*indicator of mortality*/ gen xx=1 if DementiaTenureAA-1.15<0 & (country==1 | country==3) & Dementia==0
/*indicator of mortality*/ replace xx=1 if DementiaTenureAA-4.22<0 & country==2 & Dementia==0
replace DementiaTenureAA=DementiaTenureAA-1.15 if (country==1 | country==3) & Dementia==0 & xx!=1
replace DementiaTenureAA=DementiaTenureAA-4.22 if country==2 & Dementia==0 & xx!=1
sum DementiaTenureA DementiaTenureAA
drop xx


*=================================Incident Alzheimer's disease & Vascular dementia 
drop pDementia_0 iDementia iDemAlz iDemVas
la var Dementia "All-cause dementia"
la var DemAlz "Alzheimer's disease"
la var DemVas "Vascular dementia"

gen pDementia_0=Dementia_Prevalent

tab Dementia, missing
gen iDementia=Dementia
recode iDementia 1=. if pDementia_0==1
la var iDementia "Incident all-cause dementia"
tab Dementia iDementia, missing
tab iDementia

tab DemAlz
gen iDemAlz=DemAlz
recode iDemAlz 1=. if pDementia_0==1
la var iDemAlz "Incident Alzheimer's disease"
tab DemAlz iDemAlz, missing
tab iDemAlz

tab DemVas
gen iDemVas=DemVas
recode iDemVas 1=. if pDementia_0==1
la var iDemVas "Incident vascular dementia"
tab DemVas iDemVas, missing
tab iDemVas


*=================================adding BP, Asthma, and MI
*BP
merge 1:1 IID using "data\ukb_BP0.dta"					/*53 not matched from master*/
drop _merge
la var SBPm_0 "Systolic blood pressure (mean) @baseline"
la var DBPm_0 "Diastolic blood pressure (mean) @baseline"

*Asthma&MI
merge 1:1 IID using "data\ukb_AstMI.dta"					/*53 not matched from master*/
drop _merge
gen Ast_0=1 if Ast==1 & AstBin==0 & !missing(Ast,AstBin)
recode Ast_0 .=0 if !missing(Ast,AstBin)
la var Ast_0 "Asthma @baseline"
tab Ast_0, missing

gen MI_0=1 if MI==1 & MIBin==0 & !missing(MI,MIBin)
recode MI_0 .=0 if !missing(MI,MIBin)
la var MI_0 "MI @baseline"
tab MI_0, missing

la var IHD_0 "IHD @baseline"
tab IHD_0 MI_0


*=================================SNPs
use "data\ukb_coloclised_rs6120880.dta", clear
tab rs6120880_G rs6120880_T, missing

use "data\ukb_coloclised_SNPs.dta", clear
merge 1:1 IID using "data\ukb_coloclised_rs6120880.dta", keepusing(rs6120880_G rs6120880_T)		/*All 487,409 matched*/
drop _merge
drop FID PAT MAT PHENOTYPE rs2084448_HET rs2571445_HET rs9267531_HET rs539078574_HET rs2297086_HET rs75614054_HET rs28496034_HET rs113154802_HET rs1978487_HET rs11865499_HET rs138259061_HET rs237833_HET rs11275011_HET rs548092276_HET

*Recoding for consistent effect allele with Regenie
* rs2084448_C rs9267531_G rs539078574_A rs2297086_A rs75614054_T rs28496034_G rs113154802_T rs11865499_G rs138259061_AAAAG rs11275011_TCAGTTAAAGTC rs548092276_CT rs6120880_G

*rs2084448_C 
tab rs2084448_C
recode rs2084448_C 0=2 2=0, gen(rs2084448_T)
tab rs2084448_C rs2084448_T

*rs9267531_G 
recode rs9267531_G 0=2 2=0, gen(rs9267531_A)
tab rs9267531_G rs9267531_A

*rs539078574_A 
recode rs539078574_A 0=2 2=0, gen(rs539078574_AT)
tab rs539078574_A rs539078574_AT

*rs2297086_A 
recode rs2297086_A 0=2 2=0, gen(rs2297086_G)
tab rs2297086_A rs2297086_G

*rs75614054_T 
recode rs75614054_T 0=2 2=0, gen(rs75614054_C)
tab rs75614054_T rs75614054_C

*rs28496034_G 
recode rs28496034_G 0=2 2=0, gen(rs28496034_C)
tab rs28496034_G rs28496034_C

*rs113154802_T 
recode rs113154802_T 0=2 2=0, gen(rs113154802_C)
tab rs113154802_T rs113154802_C

*rs11865499_G 
recode rs11865499_G 0=2 2=0, gen(rs11865499_A)
tab rs11865499_G rs11865499_A

*rs138259061_AAAAG 
recode rs138259061_AAAAG 0=2 2=0, gen(rs138259061_A)
tab rs138259061_AAAAG rs138259061_A

*rs11275011_TCAGTTAAAGTC 
recode rs11275011_TCAGTTAAAGTC 0=2 2=0, gen(rs11275011_T)
tab rs11275011_TCAGTTAAAGTC rs11275011_T

*rs548092276_CT
recode rs548092276_CT 0=2 2=0, gen(rs548092276_C)
tab rs548092276_CT rs548092276_C

*rs6120880_G
recode rs6120880_G 0=2 2=0, gen(rs6120880_C)
tab rs6120880_G rs6120880_C

drop rs11275011_TCAGTTAAAGTC rs113154802_T rs11865499_G rs138259061_AAAAG rs2084448_C rs2297086_A rs28496034_G rs539078574_A rs548092276_CT rs75614054_T rs9267531_G rs6120880_G

label variable rs2084448_T ""
label variable rs9267531_A ""
label variable rs539078574_AT ""
label variable rs2297086_G ""
label variable rs75614054_C ""
label variable rs28496034_C ""
label variable rs113154802_C ""
label variable rs11865499_A ""
label variable rs138259061_A ""
label variable rs11275011_T ""
label variable rs548092276_C ""
label variable rs6120880_C ""

save "data\ukb_coloclised_SNPs_m.dta", replace

use "data\ukb_LF0-DmCF_m.dta", clear
merge 1:1 IID using "data\ukb_coloclised_SNPs_m.dta"		/*15,211 not matched from master*/
drop if _merge==2 											/*156 observations deleted*/
drop _merge			


*=================================Inclusion
/*gen inc_Pairs=!missing(PairsMatching_0_logrz)
tab inc_Pairs
count if inc_Pairs==0 & !missing(ReactionTime_0_logrz)		/*395*/
count if inc_Pairs==0 & !missing(ReactionTime_2_logrz)		/*173*/
count if inc_Pairs==0 & !missing(FluidIntelligence_0_z)		/*0*/
count if inc_Pairs==0 & !missing(FluidIntelligence_2_z)		/*172*/
drop inc_Pairs */


gen inc_Reg= 1 if !missing(PairsMatching_0_logrz) | !missing(PairsMatching_0_logrz) | !missing(ReactionTime_0_logrz) | !missing(ReactionTime_2_logrz) | !missing(ProspectiveMemory_0_b) | !missing(ProspectiveMemory_2_b) | !missing(FluidIntelligence_0_z) | !missing(FluidIntelligence_2_z) | !missing(NumericMemory_0_z) | !missing(NumericMemory_2_z) | !missing(TrailMakingN_du_2_logrz) | !missing(TrailMakingAN_du_2_logrz) | !missing(TrailMakingANN_du_2_rz) | !missing(TrailMakingANNr_du_2_logrz) | !missing(SymbolDigit_c_2_z) | !missing(MatrixPattern_2_z) | !missing(TowerRearranging_c_2_z) | !missing(PairedLearning_2_z) 
recode inc_Reg .=0
count if inc_Reg==0 & !missing(ReactionTime_0_logrz)		/*0*/
count if inc_Reg==0 & !missing(ReactionTime_2_logrz)		/*0*/
count if inc_Reg==0 & !missing(FluidIntelligence_0_z)		/*0*/

la var inc_Reg "Included in Regenie analysis"
tab inc_Reg


*=================================SNPs for Sheena
keep IID age sex array centre pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 rs2571445_A rs1978487_C rs237833_A rs6120880_T rs2084448_T rs9267531_A rs539078574_AT rs2297086_G rs75614054_C rs28496034_C rs113154802_C rs11865499_A rs138259061_A rs11275011_T rs548092276_C rs6120880_C inc_Reg

*drop rs6120880_T

save "data\ukb_coloclised_SNPs_forSW.dta", replace
export delimited using "data\ukb_coloclised_SNPs_forSW", replace




*=================================Lancet commision covariates
use "data\ukb_Cov02_LF", clear
rename EID59138 IID
sort IID

save "data\ukb_Cov02_LF",replace

use "data\ukb_LF0-DmCF_m.dta", clear
merge 1:1 IID using "data\ukb_Cov02_LF"		/*108 not matched from master*/
drop _merge			

*==BMI
count if BMI_0<=15 | BMI_0>50						/*3,920 to missing*/
recode BMI_0 min/15=. 50/max=., gen(mBMI_0)			/*705 to missing*/
sum BMI_0 mBMI_0
hist mBMI_0
*BMI 3-cat
egen mBMI_0_g3=cut(mBMI_0), at(15, 25, 30, 60) icode
recode mBMI_0_g3 .=3
tabstat mBMI_0, st(n min max) by(mBMI_0_g3)
la var mBMI_0_g3 "BMI, kg/m2"
la def bmi3 0"<25" 1"25-30" 2">30" 3"Missing"
la val mBMI_0_g3 bmi3 
tab mBMI_0_g3
*BMI binary
recode mBMI_0_g 2=1 3=2, gen(mBMI_0_b)
la var mBMI_0_b "BMI, kg/m2"
la def bmi2 0"<25" 1">=25" 2"Missing"
la val mBMI_0_b bmi2 
tab mBMI_0_b
*BMI 4-cat
egen mBMI_0_g4=cut(mBMI_0), at(15, 18.5, 25, 30, 60) icode
recode mBMI_0_g4 .=4
tabstat mBMI_0, st(n min max) by(mBMI_0_g4)
la var mBMI_0_g4 "BMI, kg/m2"
la def bmi4 0"<18.5" 1"18.5-24.9" 2"25-29.9" 3">30" 4"Missing"
la val mBMI_0_g4 bmi4 
tab mBMI_0_g4

*==Birth Weight
sum Bweight_0 rBweight_0 Bweight_2 rBweight_2
drop rBweight_0 rBweight_2
compare Bweight_0 Bweight_2 if !missing(Bweight_0) & !missing(Bweight_2)
gen Bweight=Bweight_0
replace Bweight=Bweight_2 if missing(Bweight)
replace Bweight=. if !(Bweight==Bweight_2) & !missing(Bweight, Bweight_2)
la var Bweight "Birth weight"
sum Bweight, d


*===========Save
save "data\ukb_LF0-DmCF_m.dta", replace






