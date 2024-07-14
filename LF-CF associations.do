*===========================================================================
* Project: Shared Developmental Pathway (UKB)
* Purpose: Models for LF-CF associations (Lung Function-Cognitive Function)
* Launch date: 08 Apr 2022 
* Last update: 04 June 2024
*===========================================================================



/*Personal*/ cd "C:\Users\mtala\OneDrive - Queen Mary, University of London\DevOr\"
use "data\ukb_LF0-DmCF_m.dta", clear
/*QMUL*/ cd "C:\Users\hmy431\OneDrive - Queen Mary, University of London\DevOr\"
use "data\ukb_LF0-DmCF_m.dta", clear


*===Follow-up duration (FVC-Demetia)
do scripts\Models_FVC-CF.do
stset DementiaTenureAA, id(IID) failure(iDementia==1)
stcox i.FVC_b0q5 $model_4, cformat(%9.2f)				/*4,294 events; 3,904,677 PY; 320,523 subjects*/
predict HR if e(sample)==1, hr
misstable sum HR										/*320,523 subjects*/
tab iDementia if !missing(HR)							/*4,294*/
sum DementiaTenureAA if !missing(HR), d					/*mean 12.2y median 12.5 y*/
tab iDemVas if !missing(HR)								/*947*/
tab iDemAlz if !missing(HR)								/*1,855*/


*======Cases/person-years
*Note: Cases/PY by automated Tables are based on crude model 
do scripts\Models_FVC-CF.do
set more off
log using "outputs\CasesPY_Model4.log", replace
foreach oo in iDementia iDemVas iDemAlz {
	display " "
	display " "
	local ool: variable label `oo'
	display "************************************************* " "`ool'"
	display "************************************************* " 
	quietly stset DementiaTenureAA, id(IID) failure(`oo'==1)
foreach vv in FVC_b0q5 FEV1FVC_b0q5 {
	display " "
	local vvl: variable label `vv'
	display "================================================= " "`vvl'"
	quietly stcox i.`vv' $model_4, cformat(%9.2f) 
	quietly predict HR if e(sample)==1, hr
	stptime if !missing(HR), per(1000) by(`vv') dd(1)
	drop HR
}
}
log close 


*======Assumptions
*===Proportionality assumption for cox regression 
$model_4
stset DementiaTenureAA, id(IID) failure(iDementia==1)
stset DementiaTenureAA, id(IID) failure(iDemVas==1)
stset DementiaTenureAA, id(IID) failure(iDemAlz==1)
*FVC
quietly stcox i.FVC_b0q5, cformat(%9.2f)
/*Schoenfeld residuals: the null hypothesis of zero slope*/
estat phtest, detail 						/*P=0.63 0.75 0.19*/
/*Log-log plot of survival*/
stphplot, by(FVC_b0q5) plot1(msize(vtiny)) plot2(msize(vtiny)) plot3(msize(vtiny)) plot4(msize(vtiny)) plot5(msize(vtiny))

*Ratio
quietly stcox i.FEV1FVC_b0q5, cformat(%9.2f)
estat phtest, detail 						/*P=0.03 0.02 0.34*/
stphplot, by(FEV1FVC_b0q5) plot1(msize(vtiny)) plot2(msize(vtiny)) plot3(msize(vtiny)) plot4(msize(vtiny)) plot5(msize(vtiny))

*===Linear regression assumption
*FVC
set more off
log using "outputs\Cog_Fig\Assumptions\FVC_HeteroSke_MultiCol.log", replace
foreach oo in PairsMatching_0_logrz ReactionTime_0_logrz FluidIntelligence_0_z NumericMemory_0_z TrailMakingANN_du_2_rz SymbolDigit_c_2_z MatrixPattern_2_z TowerRearranging_c_2_z PairedLearning_2_z gF4_0 gF5i_2 gF9i_2 {
display "*************************************************** " "`oo'"
quietly regress `oo' i.FVC_b0q5 $model_4, cformat(%9.2f)
quietly predict r, resid
kdensity r, normal name(KD, replace) 
sleep 1000
pnorm r,  name(PP, replace)
sleep 1000
qnorm r,  name(QQ, replace)
quietly drop r
sleep 1000
rvfplot, yline(0) name(RvF, replace)
sleep 1000
graph combine KD PP QQ RvF
*quietly graph save Graph "outputs\Cog_Fig\Assumptions\FVC_`oo'.gph", replace
quietly graph export "outputs\Cog_Fig\Assumptions\FVC_`oo'.tif", replace width(1500)
estat hettest
estat vif
}
log close

*Ratio
set more off
log using "outputs\Cog_Fig\Assumptions\Ratio_HeteroSke_MultiCol.log", replace
foreach oo in PairsMatching_0_logrz ReactionTime_0_logrz FluidIntelligence_0_z NumericMemory_0_z TrailMakingANN_du_2_rz SymbolDigit_c_2_z MatrixPattern_2_z TowerRearranging_c_2_z PairedLearning_2_z gF4_0 gF5i_2 gF9i_2 {
display "*************************************************** " "`oo'"
quietly regress `oo' i.FEV1FVC_b0q5 $model_4, cformat(%9.2f)
quietly predict r, resid
kdensity r, normal name(KD, replace) 
sleep 1000
pnorm r,  name(PP, replace)
sleep 1000
qnorm r,  name(QQ, replace)
quietly drop r
sleep 1000
rvfplot, yline(0) name(RvF, replace)
sleep 1000
graph combine KD PP QQ RvF
*quietly graph save Graph "outputs\Cog_Fig\Assumptions\Ratio_`oo'.gph", replace
quietly graph export "outputs\Cog_Fig\Assumptions\Ratio_`oo'.tif", replace width(1500)
estat hettest
estat vif
}
log close

regress ReactionTime_0_logrz FVC_b0 $model_2, cformat(%9.2f)

rvpplot FVC_b0, yline(0) 	/*graphs a residual-versus-predictor plot*/
rvfplot, yline(0)			/*graphs residual-versus-fitted plot: Homoscedasticity of Residuals*/
estat hettest				/*Breusch-Pagan/Cook-Weisberg test for heteroskedasticity*/
twoway (scatter ReactionTime_0_logrz FVC_b0, msize(vtiny)) (lfit ReactionTime_0_logrz FVC_b0) (lowess ReactionTime_0_logrz FVC_b0)
							/*Linearity assumption for simple regression*/
							


*======Tables 2&3: LF-CF
do scripts\Models_FVC-CF.do
set more off
log using "outputs\LF-CF_T2&3.log", replace
foreach oo in PairsMatching_0_logrz ReactionTime_0_logrz FluidIntelligence_0_z NumericMemory_0_z TrailMakingANN_du_2_rz SymbolDigit_c_2_z MatrixPattern_2_z TowerRearranging_c_2_z PairedLearning_2_z gF4_0 gF5i_2 gF9i_2  {
	display " "
	display " "
	local ool: variable label `oo'
	display "************************************************* " "`ool'"
	display "************************************************* " 
foreach vv in FVC_b0q5 FEV1FVC_b0q5 {
	display " "
	local vvl: variable label `vv'
	display "================================================= " "`vvl'"
	regress `oo' i.`vv' $model_4, cformat(%9.2f)
forvalues i = 1(1)5  {
display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " "model `i'"
local model: display "model_" `i'
regress `oo' i.`vv' $`model', cformat(%9.2f)
display " "
}
}
}
log close 


*======Table 4: LF-Dementia
*==Cases/person-years
do scripts\Models_FVC-CF.do
set more off
log using "outputs\CasesPY_Model4.log", replace
foreach oo in iDementia iDemVas iDemAlz {
	display " "
	display " "
	local ool: variable label `oo'
	display "************************************************* " "`ool'"
	display "************************************************* " 
	quietly stset DementiaTenureAA, id(IID) failure(`oo'==1)
foreach vv in FVC_b0q5 FEV1FVC_b0q5 {
	display " "
	local vvl: variable label `vv'
	display "================================================= " "`vvl'"
	quietly stcox i.`vv' $model_4, cformat(%9.2f) 
	quietly predict HR if e(sample)==1, hr
	stptime if !missing(HR), per(1000) by(`vv') dd(1)
	drop HR
}
}
log close 

*==Association
do scripts\Models_FVC-CF.do
set more off
log using "outputs\LF-CF_T2&3.log", replace
foreach oo in iDementia iDemVas iDemAlz  {
	display " "
	display " "
	local ool: variable label `oo'
	display "************************************************* " "`ool'"
	display "************************************************* " 
foreach vv in FVC_b0q5 FEV1FVC_b0q5 {
	display " "
	local vvl: variable label `vv'
	display "================================================= " "`vvl'"
	regress `oo' i.`vv' $model_4, cformat(%9.2f)
forvalues i = 1(1)5  {
display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " "model `i'"
local model: display "model_" `i'
logistic `oo' i.`vv' $`model', cformat(%9.2f)
display " "
}
}
}
log close 


*====== Figure 2: cumulative hazard estimates after Cox regression
do scripts\Models_FVC-CF.do

*FVC
stset DementiaTenureAA, id(IID) failure(iDementia==1)
stcox i.FVC_b0q5 $model_4, cformat(%9.2f)

stcurve, cumhaz ///
	at1(FVC_b0q5=0) at2(FVC_b0q5=1) at3(FVC_b0q5=2) at4(FVC_b0q5=3) at5(FVC_b0q5=4) range(0 14.5) ///
	lcolor(cranberry orange_red navy blue midblue) lpattern(solid dash solid dash solid) ///
	yscale(range(0 0.01)) ylabel(0(0.002)0.01) ytitle(Cumulative hazard for {bf:dementia}) ylabel(, angle(horizontal) nogrid format(%9.3f)) ///
	xtitle(Years since lung function measurements) xlabel(0(2)14, nogrid) ///
	legend(order(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5") rows(5) region(lwidth(none)) title(Quintiles of {bf:FVC}, size(medsmall)) position(11) ring(0)) ///
	graphregion(fcolor(white)) title(, color(none))
graph display, ysize(5.5) xsize(7.5)
	
graph save Graph "outputs\Cog_Fig\FVCq5_Dem_CumHz.gph", replace
graph export "outputs\Cog_Fig\FVCq5_Dem_CumHz.tif", replace width(1500)


*Ratio
stcox i.FEV1FVC_b0q5 $model_4, cformat(%9.2f)

stcurve, cumhaz ///
	at1(FEV1FVC_b0q5=0) at2(FEV1FVC_b0q5=1) at3(FEV1FVC_b0q5=2) at4(FEV1FVC_b0q5=3) at5(FEV1FVC_b0q5=4) range(0 14.5) ///
	lcolor(cranberry orange_red navy blue midblue) lpattern(solid dash solid dash solid) ///
	yscale(range(0 0.01)) ylabel(0(0.002)0.01) ytitle(Cumulative hazard for {bf:dementia}) ylabel(, angle(horizontal) nogrid format(%9.3f)) ///
	xtitle(Years since lung function measurements) xlabel(0(2)14, nogrid) ///
	legend(order(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5") rows(5) region(lwidth(none)) title(Quintiles of {bf:FEV{subscript:1}/FVC}, size(medsmall)) position(11) ring(0)) ///
	graphregion(fcolor(white)) title(, color(none))
graph display, ysize(5.5) xsize(7.5)

graph save Graph "outputs\Cog_Fig\Ratioq5_Dem_CumHz.gph", replace
graph export "outputs\Cog_Fig\Ratioq5_Dem_CumHz.tif", replace width(1500)	
