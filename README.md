# sars-cov-2_genetic_drift
 Code for reproducing analysis in manuscript: 
 
 >QinQin Yu, Joao Ascensao, Takashi Okada, Olivia Boyd, Erik Volz, and Oskar Hallatschek, Lineage frequency time series reveal elevated levels of genetic drift in SARS-CoV-2 transmission in England, 2022.
 
# Abstract 

# Description of repository

 In this manuscript, we develop and validate a method for jointly inferring time-varying effective population size and measurement noise overdispersion from time series lineage frequency data. The model using a Hidden Markov Model (HMM) where the hidden states are the true frequencies over time, and the observed states are the observed frequencies. The easiest way to run this code is to download the repository and run
 
 `python analysis_scripts/run_hmm_inference/run_hmm_inference.py`
 
 The rest of the repository contains all of the code and data needed for reproducing the additional analyses described in the paper. 
 
### analysis_scripts/
This folder contains all analysis scripts used for analyzing the actual data of lineage frequency time series in England. 
- cut_tree/: Code for cutting the phylogenetic tree to create lineages.
- get_lineage_counts/: Code for getting the number of sequences in each lineage (data reformatting). 
- run_hmm_inference/: Code for running the inference of effective population size and measurement noise overdispersion. 
- fitness/: Code for jointly inferring the selection coefficient of each lineage and the population-wide effective population size. 
- epidemiological_models/: Code for calculating the effective population size for SIR and SEIR compartmental models using data of the number of infected individuals and the effective reproduction number.

### functions/
Functions called by analysis and simulation codes. 

### data/
Data from England used in analysis. 

### simulations/
Code and output of 3 types of simulations: 
- simulation_scripts/wf_simulation.py: Wright-Fisher simulations with observation step for validating the method
- simulation_scripts/stochastic_seir_simulation/: Stochastic SEIR model for validating calculation of SEIR model effective population size
- simulation_scripts/deme_simulations/: Simulations of deme structure for testing its effect on the effectie population size. 

### figures/
Code for reproducing all figures in manuscript and figure outputs.

### manuscript/
Miscellaneous files used for manuscript submission. 


 
 
