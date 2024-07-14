*===========================================================================
* Project: Shared Developmental Pathway (UKB)
* Purpose: Models for LuFu-CogFu associations
* Launch date: 08 Apr 2022 
* Last date: 08 Apr 2022
* Changes: 
*===========================================================================

*Defining the model macros:
global model_1 " "
global model_2 "$model_1 age i.sex i.centre_m TDI"
global model_3 "$model_2 i.Edu_0m i.LC2_Hearing_0 i.LC3_HeadInjury_0 i.LC4_HTNmLife_0 i.LC5_Alcohol_0 i.LC6_Obese_0 i.LC7_T2D_0 i.Smoking_0m i.LC9_SocialIso_0 i.LC10_Depression_0 i.LC11_PhInactivity_0 i.LC12_AirPollution"
global model_4 "$model_3 i.APOE4_alleles"
global model_5 "$model_4 Bweight" 
*===========================================================================

