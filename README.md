# MCMV-primary-infection
Found here is the R code used to analyze data presented in the article titled "Spatial kinetics and immune control of murine cytomegalovirus infection in the salivary glands". All code was written by Catherine Byrne. All queries can be directed to cbyrne@fredhutch.org.

If wanting to reproduce results in the article, all data and R files should first be downloaded. 

Data files are described below.
1. mcmv_sg_inoc_flow_results.csv - contains all results from staining and performing flow cytometry on the blood, spleens, and salivary glands of mice infected with MCMV via the salivary glands.
2. mcmv_sg_inoc_luminescence_results.csv - contains results from daily bioimaging of mice that were infected with MCMV via the salivary glands.
3. mcmv_ip_inoc_flow_results.csv - contains resuts from staining and performing flow cytometry on the blood of mice infected with MCMV via an intraperitoneal infection.
4. mcmv_ip_inoc_luminescence_results.csv - contains results from daily bioimaging of mice that were infected with MCMV via an intraperitoneal injection.
   
R files can be run in the following order
1. combined_cytokines_immune_plot.r - This file compares how well a model with both cytokines and virus-specific CD8s vs just virus-specific CD8s clearing MCMV infection fits data. 
2. mcmv_pomp_cytokines_analyze.r - This file allows the model with both cytokines and virus-specific CD8s to be fit to all mouse data and summarizes results.
3. transient_infection_analysis.r - This file analyzes the predicted viral dynamics when low doses of virus are given to mice. It predicts under what conditions a sucessful, unsuccessful, or transient infection will occur.

Please note, the data files mcmv_sg_inoc_flow_results.csv and mcmv_sg_inoc_luminescence_results.csv  were used to produce the main results in the paper. The code can be edited to study how results are predicted to change when viral inoculations are delivered through intraperitoneal injection. Further, because this code fits data using a stochastic method, results may slightly differ from those found in the paper.
