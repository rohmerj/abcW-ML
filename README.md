# abcW-ML

This repository contains:

- [ ] R scripts
  run_doe.R : generate a deisgn of experiments from cyclone database
  run_abcW_cyclone.R : run the analysis for varying proportions of training samples (with repetitions)
  run_valid_ML.R : run the ML validation procedure

- [ ] Files (please unzip them before launching the scripts)
  ./data: Cyclone characteristics and Hs extracted at the two points
  ./doe: design of experiments for varying proportions of training samples (with repetitions)
  ./results: EVA results for varying proportions and different sampling procedures. It contains the R script for reproducing Figs. 6-9
  ./valid: ML validation results for varying proportions. It contains the R script for reproducing Fig. 5 and Suppl. Mat. B 

The necessary R libraries are:
ranger
POT
transport
clhs
ggplot2
gridExtra
dplyr
