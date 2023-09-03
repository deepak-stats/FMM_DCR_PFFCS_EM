# FMM_DCR_PFFCS_EM

This folder contain the dataset and simulation codes:

DataSet codes:

Cause_2_EP_data_EM.R 

Cause_2_GE_data_EM.R

Cause_3_EP_data_EM.R

Cause_3_GE_data_EM.R


For two causes the  dataset is given in Cause_2.csv and for three causes the dataset is given in Cause_3.csv. 

Simulation codes:

EP_distribution_EM.R

GE_distribution_EM.R

EP_CPs and ALs.R

GE_CPs and ALs.R

Setting the paramter values, n, m, k and the scheme will give you the average estimates, mean square errors, coverage probabilitys (CPs) and average lenghths (ALs). For example, the settings are set in the code as: alpha_1=4.75, lambda_1=2.25, alpha_2=1.75, lambda_2=1.20, pi=0.40, k=3, n=60, m=60, with scheme III R=c(rep(0,m/2-1),(n-m),rep(0,m/2)). Feel free to contact me if you have any problem in running the codes, Email id: deepakprajapati@iiml.ac.in. 

