# MCMV-primary-infection
Found here is the R code used to analyze data presented in the article titled "Spatial kinetics and immune control of murine cytomegalovirus infection in the salivary glands". All code was written by Catherine Byrne. All queries can be directed to cbyrne@fredhutch.org.

If wanting to reproduce results in the article, all data and R files should first be downloaded. 

Data files are described below.
1. all_flow_results_new.csv - contains all results from staining and performing flow cytometry on the blood, spleens, and salivary glands of mice infecte with MCMV via the salivary glands.
2. all_luminescence_results_new.csv - contains results from daily bioimaging of mice.
   
R files can be run in the following order
1. combined_cytokines_immune_plot.r - This file compares how well a model with both cytokines and virus-specific CD8s vs just virus-specific CD8s clearing MCMV infection fits data. 
2. mcmv_pomp_cytokines_analyze.r - This file allows the model with both cytokines and virus-specific CD8s to be fit to all mouse data and summarizes results.
3. transient_infection_analysis.r - This file analyzes the predicted viral dynamics when low doses of virus are given to mice. It predicts under what conditions a sucessful, unsuccessful, or transient infection will occur.
