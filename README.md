# sars-cov-2_genetic_drift
 Code for reproducing analysis in manuscript: 
 
**Title**: 
Lineage frequency time series reveal elevated levels of genetic drift in SARS-CoV-2 transmission in England [[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.11.21.517390v1)]

**Authors**: 
QinQin Yu, Joao Ascensao, Takashi Okada, Olivia Boyd, Erik Volz, and Oskar Hallatschek
 
**Abstract**:
Random genetic drift in the population-level dynamics of an infectious disease outbreak results from the randomness of inter-host transmission and the randomness of host recovery or death. The strength of genetic drift has been found to be high for SARS-CoV-2 due to superspreading, and this is expected to substantially impact the disease epidemiology and evolution. Noise that results from the measurement process, such as biases in data collection across time, geographical areas, etc., can potentially confound estimates of genetic drift as both processes contribute "noise" to the data. To address this challenge, we develop and validate a method to jointly infer genetic drift and measurement noise from time-series lineage frequency data. We apply this method to over 490,000 SARS-CoV-2 genomic sequences from England collected between March 2020 and December 2021 by the COVID-19 Genomics UK (COG-UK) consortium. We find that even after correcting for measurement noise, the strength of genetic drift is consistently, throughout time, higher than that expected from the observed number of COVID-19 positive individuals in England by 1 to 3 orders of magnitude. Corrections taking into account epidemiological dynamics (susceptible-infected-recovered or susceptible-exposed-infected-recovered models) do not explain the discrepancy. Moreover, the levels of genetic drift that we observe are higher than the estimated levels of superspreading found by modeling studies that incorporate data on actual contact statistics in England. We discuss how even in the absence of superspreading, high levels of genetic drift can be generated via community structure in the host contact network. Our results suggest that further investigations of heterogeneous host contact structure may be important for understanding the high levels of genetic drift observed for SARS-CoV-2 in England.

# Description of repository

**Note that this repository contains some large files (>100MB, specifically the tree files) that require [Git Large File Storage](https://git-lfs.com/) to be installed in order to clone.**

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


 
 
