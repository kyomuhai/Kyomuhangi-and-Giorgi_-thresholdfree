## Kyomuhangi-and-Giorgi_-thresholdfree
## Title:  A threshold-free approach with age-dependency for estimating malaria seroprevalence

authors:  Irene Kyomuhangi (i.kyomuhangi@lancaster.ac.uk, kyomuhai@gmail.com)  and Emanuele Giorgi, Lancaster University (e.giorgi@lancaster.ac.uk), CHICAS,  Lancaster University 


# Background
This script contains syntax used for the analysis of malaria serology data as described in the paper: 
"A threshold-free approach with age-dependency for estimating malaria seroprevalence", as specified for the M2 approach. 
For simplicity, code for M1 is not provided, however this is available on request from the authors. 
The antibody measurements used in this analysis are PfAMA OD values obtained from ELISA, however the methods are applicable to any malaria antigen type, and continuous antibody measurement. 

The code is split into 2 parts:
1) Part A implements M2, the threshold-free approach introduced in this paper, and 
2) Part B implements M1, the classic threshold dependent approach as described in the paper. 

Throughout the script we indicate where the code may need to change depending on the dataset/antibody type under analysis. 
Explanations of the functions and operations used are provided within the syntax itself, and further details of statistical/mathematical principles of this analysis can be found in the paper. 

To request access to the dataset used, please contact Gillian Stresman (Gillian.Stresman@lshtm.ac.uk) or Chris Drakeley (Chris.Drakeley@lshtm.ac.uk) at the London School of Hygiene and Tropical Medicine.  


# Packages
This code is run in R. Packages to install: numDeriv, tidyverse, latex2exp, data.table, lme4, stringr, dplyr, haven

# Data

Data for this analysis should contain:
1) continuous antibody measurements (eg, OD and MFI). Note that this analysis does not use seropositive/seronegative values
2) age in years

Ensure that your dataset does not contain missing values for either of these variables. 

This code can be extended to include analysis using other variables/covariates 9e.g altitude), and we have indicated where this is possible in the syntax. 
