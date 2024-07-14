*===========================================================================
* Project: Shared Developmental Pathway (UKB) 
* Purpose: Table 1
* Output files: 
* Launch date: 02 March 2023 
*===========================================================================

/*QMUL*/


/*Home*/ cd "C:\Users\mtala\OneDrive - Queen Mary, University of London\DevOr\"
use "data\ukb_LF0-DmCF_m.dta", clear
/*Office*/ cd "C:\Users\hmy431\OneDrive - Queen Mary, University of London\DevOr\"
use "data\ukb_LF0-DmCF_m.dta", clear


cd "C:\Users\mtala\Documents" 


*============================== FVC header =================================*
*count if AstpLuFu==1 & !missing(kq930q4) & kq930q4 !=4
putexcel set outputs\Table1_DevOr.xlsx, sheet(T1_FVC) modify

*-->Defining culumn variable
*gen colvar=kq930q4 if AstpLuFu==1 & kq930q4 !=4
gen colvar=FVC_b0q5
putexcel C1=("Quintiles of FVC"), bold

putexcel B2=("Q1") C2=("Q2") D2=("Q3") E2=("Q4") F2=("Q5") G2=("Q6"), bold

*n
putexcel A3=("n"), bold
tabulate colvar, matcell(nn)
mat l nn
putexcel B3 = nn[1,1] C3 = nn[2,1] D3 = nn[3,1] E3 = nn[4,1] F3 = nn[5,1] G3 = nn[6,1], names nformat("#,###")

*Continuous 
gen convar=FVC_b0
putexcel A4= ("FVC (L) best measure"), bold

tabstat convar, stat(mean sd) by(colvar) col(stat) long format(%9.3g) save
matrix Q1 = r(Stat1) 
matrix Q2 = r(Stat2)
matrix Q3 = r(Stat3)
matrix Q4 = r(Stat4)
matrix Q5 = r(Stat5)
matrix Q6 = r(Stat6)
local msd1 : display %3.2f Q1[1,1] " ± " %3.1f Q1[2,1]
local msd2 : display %3.2f Q2[1,1] " ± " %3.1f Q2[2,1]
local msd3 : display %3.2f Q3[1,1] " ± " %3.1f Q3[2,1]
local msd4 : display %3.2f Q4[1,1] " ± " %3.1f Q4[2,1]
local msd5 : display %3.2f Q5[1,1] " ± " %3.1f Q5[2,1]
local msd6 : display %3.2f Q6[1,1] " ± " %3.1f Q6[2,1]
putexcel B4 = ("`msd1'") C4 = ("`msd2'") D4 = ("`msd3'") E4 = ("`msd4'") F4 = ("`msd5'") G4 = ("`msd6'"), names nformat(number_d2) 
drop convar


*============================== Ratio header =================================*
putexcel set outputs\Table1_DevOr.xlsx, sheet(T1_Ratio) modify

*-->Defining culumn variable
gen colvar=FEV1FVC_b0q5
putexcel C1=("Quintiles of FEV1/FVC ratio"), bold

putexcel B2=("Q1") C2=("Q2") D2=("Q3") E2=("Q4") F2=("Q5") G2=("Q6"), bold

*n
putexcel A3=("n"), bold
tabulate colvar, matcell(nn)
mat l nn
putexcel B3 = nn[1,1] C3 = nn[2,1] D3 = nn[3,1] E3 = nn[4,1] F3 = nn[5,1] G3 = nn[6,1], names nformat("#,###")

*Continuous 
gen convar=FEV1FVC_b0
putexcel A4= ("FEV1/FVC ratio best measure"), bold

tabstat convar, stat(mean sd) by(colvar) col(stat) long format(%9.3g) save
matrix Q1 = r(Stat1) 
matrix Q2 = r(Stat2)
matrix Q3 = r(Stat3)
matrix Q4 = r(Stat4)
matrix Q5 = r(Stat5)
matrix Q6 = r(Stat6)
local msd1 : display %3.2f Q1[1,1] " ± " %3.1f Q1[2,1]
local msd2 : display %3.2f Q2[1,1] " ± " %3.1f Q2[2,1]
local msd3 : display %3.2f Q3[1,1] " ± " %3.1f Q3[2,1]
local msd4 : display %3.2f Q4[1,1] " ± " %3.1f Q4[2,1]
local msd5 : display %3.2f Q5[1,1] " ± " %3.1f Q5[2,1]
local msd6 : display %3.2f Q6[1,1] " ± " %3.1f Q6[2,1]
putexcel B4 = ("`msd1'") C4 = ("`msd2'") D4 = ("`msd3'") E4 = ("`msd4'") F4 = ("`msd5'") G4 = ("`msd6'"), names nformat(number_d2) 
drop convar


*=======================================Body of table
*=======================================
*-->Age
gen convar=age
putexcel A5= ("Age, y"), bold

tabstat convar, stat(mean sd) by(colvar) col(stat) long format(%9.3g) save
matrix Q1 = r(Stat1) 
matrix Q2 = r(Stat2)
matrix Q3 = r(Stat3)
matrix Q4 = r(Stat4)
matrix Q5 = r(Stat5)
matrix Q6 = r(Stat6)
local msd1 : display %3.1f Q1[1,1] " ± " %3.1f Q1[2,1]
local msd2 : display %3.1f Q2[1,1] " ± " %3.1f Q2[2,1]
local msd3 : display %3.1f Q3[1,1] " ± " %3.1f Q3[2,1]
local msd4 : display %3.1f Q4[1,1] " ± " %3.1f Q4[2,1]
local msd5 : display %3.1f Q5[1,1] " ± " %3.1f Q5[2,1]
local msd6 : display %3.1f Q6[1,1] " ± " %3.1f Q6[2,1]
putexcel B5 = ("`msd1'") C5 = ("`msd2'") D5 = ("`msd3'") E5 = ("`msd4'") F5 = ("`msd5'") G5 = ("`msd6'"), names nformat(number_d2) 
drop convar


*-->Sex
gen binvar=sex
putexcel A6= ("Men, n (%)"), bold

tabulate binvar colvar, col matcell(nn) 
local np1 : display %3.0f nn[2,1] " (" %3.1f (nn[2,1]/(nn[2,1]+nn[1,1]))*100 ")"
local np2 : display %3.0f nn[2,2] " (" %3.1f (nn[2,2]/(nn[2,2]+nn[1,2]))*100 ")"
local np3 : display %3.0f nn[2,3] " (" %3.1f (nn[2,3]/(nn[2,3]+nn[1,3]))*100 ")"
local np4 : display %3.0f nn[2,4] " (" %3.1f (nn[2,4]/(nn[2,4]+nn[1,4]))*100 ")"
local np5 : display %3.0f nn[2,5] " (" %3.1f (nn[2,5]/(nn[2,5]+nn[1,5]))*100 ")"
local np6 : display %3.0f nn[2,6] " (" %3.1f (nn[2,6]/(nn[2,6]+nn[1,6]))*100 ")"
putexcel B6 = ("`np1'") C6 = ("`np2'") D6 = ("`np3'") E6 = ("`np4'") F6 = ("`np5'") G6 = ("`np6'"), names nformat(number_d2) 
drop binvar


*-->TDI
gen convar=TDI
putexcel A7= ("Townsend deprivation index"), bold

tabstat convar, stat(mean sd) by(colvar) col(stat) long format(%9.3g) save
matrix Q1 = r(Stat1) 
matrix Q2 = r(Stat2)
matrix Q3 = r(Stat3)
matrix Q4 = r(Stat4)
matrix Q5 = r(Stat5)
matrix Q6 = r(Stat6)
local msd1 : display %3.2f Q1[1,1] " ± " %3.1f Q1[2,1]
local msd2 : display %3.2f Q2[1,1] " ± " %3.1f Q2[2,1]
local msd3 : display %3.2f Q3[1,1] " ± " %3.1f Q3[2,1]
local msd4 : display %3.2f Q4[1,1] " ± " %3.1f Q4[2,1]
local msd5 : display %3.2f Q5[1,1] " ± " %3.1f Q5[2,1]
local msd6 : display %3.2f Q6[1,1] " ± " %3.1f Q6[2,1]
putexcel B7 = ("`msd1'") C7 = ("`msd2'") D7 = ("`msd3'") E7 = ("`msd4'") F7 = ("`msd5'") G7 = ("`msd6'"), names nformat(number_d2) 
drop convar


*-->Height
gen convar=Height_0
putexcel A8= ("Height (cm)"), bold

tabstat convar, stat(mean sd) by(colvar) col(stat) long format(%9.3g) save
matrix Q1 = r(Stat1) 
matrix Q2 = r(Stat2)
matrix Q3 = r(Stat3)
matrix Q4 = r(Stat4)
matrix Q5 = r(Stat5)
matrix Q6 = r(Stat6)
local msd1 : display %3.0f Q1[1,1] " ± " %3.1f Q1[2,1]
local msd2 : display %3.0f Q2[1,1] " ± " %3.1f Q2[2,1]
local msd3 : display %3.0f Q3[1,1] " ± " %3.1f Q3[2,1]
local msd4 : display %3.0f Q4[1,1] " ± " %3.1f Q4[2,1]
local msd5 : display %3.0f Q5[1,1] " ± " %3.1f Q5[2,1]
local msd6 : display %3.0f Q6[1,1] " ± " %3.1f Q6[2,1]
putexcel B8 = ("`msd1'") C8 = ("`msd2'") D8 = ("`msd3'") E8 = ("`msd4'") F8 = ("`msd5'") G8 = ("`msd6'"), names nformat(number_d2) 
drop convar


*-->Categorical with 5 levels
gen cat5var=Edu_0m5
putexcel A9= ("Qualifications, n (%)"), bold
putexcel A10= ("College or University degree"), bold txtin(2)
putexcel A11= ("A/AS levels"), bold txtin(2)
putexcel A12= ("O levels/GCSEs"), bold txtin(2)
putexcel A13= ("CSEs/NVQ/HND/HNC or Other"), bold txtin(2)
putexcel A14= ("None of the above"), bold txtin(2)

tabulate cat5var colvar, col matcell(nn) 
local np1 : display %3.0f nn[1,1] " (" %3.1f (nn[1,1]/(nn[1,1]+nn[2,1]+nn[3,1]+nn[4,1]+nn[5,1]))*100 ")"
local np2 : display %3.0f nn[1,2] " (" %3.1f (nn[1,2]/(nn[1,2]+nn[2,2]+nn[3,2]+nn[4,2]+nn[5,2]))*100 ")"
local np3 : display %3.0f nn[1,3] " (" %3.1f (nn[1,3]/(nn[1,3]+nn[2,3]+nn[3,3]+nn[4,3]+nn[5,3]))*100 ")"
local np4 : display %3.0f nn[1,4] " (" %3.1f (nn[1,4]/(nn[1,4]+nn[2,4]+nn[3,4]+nn[4,4]+nn[5,4]))*100 ")"
local np5 : display %3.0f nn[1,5] " (" %3.1f (nn[1,5]/(nn[1,5]+nn[2,5]+nn[3,5]+nn[4,5]+nn[5,5]))*100 ")"
local np6 : display %3.0f nn[1,6] " (" %3.1f (nn[1,6]/(nn[1,6]+nn[2,6]+nn[3,6]+nn[4,6]+nn[5,6]))*100 ")"
putexcel B10 = ("`np1'") C10 = ("`np2'") D10 = ("`np3'") E10 = ("`np4'") F10 = ("`np5'") G10 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[2,1] " (" %3.1f (nn[2,1]/(nn[1,1]+nn[2,1]+nn[3,1]+nn[4,1]+nn[5,1]))*100 ")"
local np2 : display %3.0f nn[2,2] " (" %3.1f (nn[2,2]/(nn[1,2]+nn[2,2]+nn[3,2]+nn[4,2]+nn[5,2]))*100 ")"
local np3 : display %3.0f nn[2,3] " (" %3.1f (nn[2,3]/(nn[1,3]+nn[2,3]+nn[3,3]+nn[4,3]+nn[5,3]))*100 ")"
local np4 : display %3.0f nn[2,4] " (" %3.1f (nn[2,4]/(nn[1,4]+nn[2,4]+nn[3,4]+nn[4,4]+nn[5,4]))*100 ")"
local np5 : display %3.0f nn[2,5] " (" %3.1f (nn[2,5]/(nn[1,4]+nn[2,5]+nn[3,5]+nn[4,5]+nn[5,5]))*100 ")"
local np6 : display %3.0f nn[2,6] " (" %3.1f (nn[2,6]/(nn[1,4]+nn[2,6]+nn[3,6]+nn[4,6]+nn[5,6]))*100 ")"
putexcel B11 = ("`np1'") C11 = ("`np2'") D11 = ("`np3'") E11 = ("`np4'") F11 = ("`np5'") G11 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[3,1] " (" %3.1f (nn[3,1]/(nn[1,1]+nn[2,1]+nn[3,1]+nn[4,1]+nn[5,1]))*100 ")"
local np2 : display %3.0f nn[3,2] " (" %3.1f (nn[3,2]/(nn[1,2]+nn[2,2]+nn[3,2]+nn[4,2]+nn[5,2]))*100 ")"
local np3 : display %3.0f nn[3,3] " (" %3.1f (nn[3,3]/(nn[1,3]+nn[2,3]+nn[3,3]+nn[4,3]+nn[5,3]))*100 ")"
local np4 : display %3.0f nn[3,4] " (" %3.1f (nn[3,4]/(nn[1,4]+nn[2,4]+nn[3,4]+nn[4,4]+nn[5,4]))*100 ")"
local np5 : display %3.0f nn[3,5] " (" %3.1f (nn[3,5]/(nn[1,5]+nn[2,5]+nn[3,5]+nn[4,5]+nn[5,5]))*100 ")"
local np6 : display %3.0f nn[3,6] " (" %3.1f (nn[3,6]/(nn[1,6]+nn[2,6]+nn[3,6]+nn[4,6]+nn[5,6]))*100 ")"
putexcel B12 = ("`np1'") C12 = ("`np2'") D12 = ("`np3'") E12 = ("`np4'") F12 = ("`np5'") G12 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[4,1] " (" %3.1f (nn[4,1]/(nn[1,1]+nn[2,1]+nn[3,1]+nn[4,1]+nn[5,1]))*100 ")"
local np2 : display %3.0f nn[4,2] " (" %3.1f (nn[4,2]/(nn[1,2]+nn[2,2]+nn[3,2]+nn[4,2]+nn[5,2]))*100 ")"
local np3 : display %3.0f nn[4,3] " (" %3.1f (nn[4,3]/(nn[1,3]+nn[2,3]+nn[3,3]+nn[4,3]+nn[5,3]))*100 ")"
local np4 : display %3.0f nn[4,4] " (" %3.1f (nn[4,4]/(nn[1,4]+nn[2,4]+nn[3,4]+nn[4,4]+nn[5,4]))*100 ")"
local np5 : display %3.0f nn[4,5] " (" %3.1f (nn[4,5]/(nn[1,5]+nn[2,5]+nn[3,5]+nn[4,5]+nn[5,5]))*100 ")"
local np6 : display %3.0f nn[4,6] " (" %3.1f (nn[4,6]/(nn[1,6]+nn[2,6]+nn[3,6]+nn[4,6]+nn[5,6]))*100 ")"
putexcel B13 = ("`np1'") C13 = ("`np2'") D13 = ("`np3'") E13 = ("`np4'") F13 = ("`np5'") G13 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[5,1] " (" %3.1f (nn[5,1]/(nn[1,1]+nn[2,1]+nn[3,1]+nn[4,1]+nn[5,1]))*100 ")"
local np2 : display %3.0f nn[5,2] " (" %3.1f (nn[5,2]/(nn[1,2]+nn[2,2]+nn[3,2]+nn[4,2]+nn[5,2]))*100 ")"
local np3 : display %3.0f nn[5,3] " (" %3.1f (nn[5,3]/(nn[1,3]+nn[2,3]+nn[3,3]+nn[4,3]+nn[5,3]))*100 ")"
local np4 : display %3.0f nn[5,4] " (" %3.1f (nn[5,4]/(nn[1,4]+nn[2,4]+nn[3,4]+nn[4,4]+nn[5,4]))*100 ")"
local np5 : display %3.0f nn[5,5] " (" %3.1f (nn[5,5]/(nn[1,5]+nn[2,5]+nn[3,5]+nn[4,5]+nn[5,5]))*100 ")"
local np6 : display %3.0f nn[5,6] " (" %3.1f (nn[5,6]/(nn[1,6]+nn[2,6]+nn[3,6]+nn[4,6]+nn[5,6]))*100 ")"
putexcel B14 = ("`np1'") C14 = ("`np2'") D14 = ("`np3'") E14 = ("`np4'") F14 = ("`np5'") G14 = ("`np6'"), names nformat(number_d2) 
drop cat5var


*-->Smoking
gen cat3var=Smoking_0m
putexcel A15= ("Smoking, n (%)"), bold
putexcel A16= ("Never"), bold txtin(2)
putexcel A17= ("Previous"), bold txtin(2)
putexcel A18= ("Current"), bold txtin(2)

tabulate cat3var colvar, col matcell(nn) 
local np1 : display %3.0f nn[1,1] " (" %3.1f (nn[1,1]/(nn[1,1]+nn[2,1]+nn[3,1]))*100 ")"
local np2 : display %3.0f nn[1,2] " (" %3.1f (nn[1,2]/(nn[1,2]+nn[2,2]+nn[3,2]))*100 ")"
local np3 : display %3.0f nn[1,3] " (" %3.1f (nn[1,3]/(nn[1,3]+nn[2,3]+nn[3,3]))*100 ")"
local np4 : display %3.0f nn[1,4] " (" %3.1f (nn[1,4]/(nn[1,4]+nn[2,4]+nn[3,4]))*100 ")"
local np5 : display %3.0f nn[1,5] " (" %3.1f (nn[1,5]/(nn[1,5]+nn[2,5]+nn[3,5]))*100 ")"
local np6 : display %3.0f nn[1,6] " (" %3.1f (nn[1,6]/(nn[1,6]+nn[2,6]+nn[3,6]))*100 ")"
putexcel B16 = ("`np1'") C16 = ("`np2'") D16 = ("`np3'") E16 = ("`np4'") F16 = ("`np5'") G16 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[2,1] " (" %3.1f (nn[2,1]/(nn[1,1]+nn[2,1]+nn[3,1]))*100 ")"
local np2 : display %3.0f nn[2,2] " (" %3.1f (nn[2,2]/(nn[1,2]+nn[2,2]+nn[3,2]))*100 ")"
local np3 : display %3.0f nn[2,3] " (" %3.1f (nn[2,3]/(nn[1,3]+nn[2,3]+nn[3,3]))*100 ")"
local np4 : display %3.0f nn[2,4] " (" %3.1f (nn[2,4]/(nn[1,4]+nn[2,4]+nn[3,4]))*100 ")"
local np5 : display %3.0f nn[2,5] " (" %3.1f (nn[2,5]/(nn[1,4]+nn[2,5]+nn[3,5]))*100 ")"
local np6 : display %3.0f nn[2,6] " (" %3.1f (nn[2,6]/(nn[1,4]+nn[2,6]+nn[3,6]))*100 ")"
putexcel B17 = ("`np1'") C17 = ("`np2'") D17 = ("`np3'") E17 = ("`np4'") F17 = ("`np5'") G17 = ("`np6'"), names nformat(number_d2) 
local np1 : display %3.0f nn[3,1] " (" %3.1f (nn[3,1]/(nn[1,1]+nn[2,1]+nn[3,1]))*100 ")"
local np2 : display %3.0f nn[3,2] " (" %3.1f (nn[3,2]/(nn[1,2]+nn[2,2]+nn[3,2]))*100 ")"
local np3 : display %3.0f nn[3,3] " (" %3.1f (nn[3,3]/(nn[1,3]+nn[2,3]+nn[3,3]))*100 ")"
local np4 : display %3.0f nn[3,4] " (" %3.1f (nn[3,4]/(nn[1,4]+nn[2,4]+nn[3,4]))*100 ")"
local np5 : display %3.0f nn[3,5] " (" %3.1f (nn[3,5]/(nn[1,5]+nn[2,5]+nn[3,5]))*100 ")"
local np6 : display %3.0f nn[3,6] " (" %3.1f (nn[3,6]/(nn[1,6]+nn[2,6]+nn[3,6]))*100 ")"
putexcel B18 = ("`np1'") C18 = ("`np2'") D18 = ("`np3'") E18 = ("`np4'") F18 = ("`np5'") G18 = ("`np6'"), names nformat(number_d2) 
drop cat3var

drop colvar
