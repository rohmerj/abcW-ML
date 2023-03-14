# abcW-ML: combining ABC, Wasserstein distance and ML predictions for extreme value analysis

This repository contains the data and R scripts to reproduce the results of Rohmer et al. (2023) "Combining machine learning predictions and numerical simulation results for the extreme value analysis of cyclone-induced wave heights – Application in Guadeloupe", submitted to Ocean Modelling

- [ ] R scripts
  - [run_doe.R](./run_doe.R) : generate a design of experiments from the cyclone database
  - run_abcW_cyclone.R : run the analysis for varying proportions of training samples (with repetitions)
  - run_valid_ML.R : run the ML validation procedure

- [ ] Files (please unzip them before launching the scripts)
  - ./data: Cyclone characteristics and Hs extracted at the two points
  - ./doe: design of experiments for varying proportions of training samples (with repetitions)
  - ./results: EVA results for varying proportions and different sampling procedures. It contains the R script for reproducing Figs. 6-9
  - ./valid: ML validation results for varying proportions. It contains the R script for reproducing Fig. 5 and Suppl. Mat. B 

- [ ] The necessary R libraries are:
  - [ranger](https://cran.r-project.org/web/packages/ranger/index.html)
  - [POT](https://cran.r-project.org/web/packages/POT/index.html)
  - [transport](https://cran.r-project.org/web/packages/transport/index.html)
  - [clhs](https://cran.r-project.org/web/packages/clhs/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
  - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
