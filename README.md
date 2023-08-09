# abcW-ML: combining Wasserstein distance, ABC algorithm, and uncertain ML predictions for extreme value analysis

This repository contains the data and R scripts to reproduce the results of Rohmer et al. (2023) "Combining uncertain machine learning predictions and numerical simulation results for the extreme value analysis of cyclone-induced wave heights – Application in Guadeloupe", submitted to Ocean Modelling

- [ ] R scripts
  - [run_doe_MLtraining.R](./run_doe_MLtraining.R) : generate a design of experiments from the cyclone database and run the ML training + validation procedure
  - [run_abcW_cyclone.R](./run_abcW_cyclone.R) : run the analysis for varying proportions of training samples (with repetitions)

- [ ] Files (please unzip them before launching the scripts)
  - ./data: Cyclone characteristics and Hs extracted at the four points
  - ./doe: design of experiments for varying proportions of training samples (with repetitions)
  - ./results: EVA results for varying proportions and different sampling procedures. It contains the R script for reproducing Figs. 6,7, 9 and 10
  - ./valid: ML validation results for varying proportions. It contains the R script for reproducing Fig. 5 and 8 and Suppl. Mat. C 

- [ ] The necessary R libraries are:
  - [ranger](https://cran.r-project.org/web/packages/ranger/index.html)
  - [POT](https://cran.r-project.org/web/packages/POT/index.html)
  - [transport](https://cran.r-project.org/web/packages/transport/index.html)
  - [clhs](https://cran.r-project.org/web/packages/clhs/index.html)
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
  - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
  - [scoringRules](https://cran.r-project.org/web/packages/scoringRules/index.html)
