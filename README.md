# A Robustness Test for Estimating Total Effects with Covariate Adjustment

This repository contains the code in the paper "A Robustness Test for Estimating Total Effects with Covariate Adjustment".

RoCA is an acronym for **Ro**bust **C**ovariate **A**djustment.

## File hierarchy

```
.
├── RoCA                                  # implementation of the testing procedure
│   ├── RoCA.R                            # main function of the test
│   ├── utils.R                           # utility functions for RoCA.R
│   ├── utils-vignette.R                  # utilities for vignette.R
│   └── vignette.R                        # examples for the usage of RoCA 
├── data                                  # data from Sachs et al. (2005)
│   ├── 10. cd3cd28icam2+aktinhib.xls
│   └── ...
├── simulation                            # code for the simulation
│   ├── example.R                         # code for Example 9
│   ├── R                                 # folder of R functions used in the simulation
│   │   ├── get_true_covmat.R
│   │   └── ...
│   ├── runEx.sh                          # script for running example.R
│   ├── runSim.sh                         # script for running simulation.R
│   ├── run_simulation.R                  # code for the simulation study
│   ├── sc_graph.R                        # code for the real data example
│   └── simulation.R                      # wrapper of run_simulation.R
└── README.md
```
