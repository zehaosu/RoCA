# A Robustness Test for Estimating Total Effects with Covariate Adjustment

This repository contains the code in the paper "A Robustness Test for Estimating Total Effects with Covariate Adjustment" by Zehao Su and Leonard Henckel.

RoCA is an acronym for **Ro**bust **C**ovariate **A**djustment.

## File hierarchy

```
.
├── data                                  # data from Sachs et al. (2005)
│   ├── 10. cd3cd28icam2+aktinhib.xls
│   └── ...
├── example.R                             # code for Example 9
├── R                                     # R code for the implementation of our test
│   ├── get_true_covmat.R
│   ├── get_unique_node_vas.R
│   ├── irm.R
│   ├── main.R
│   ├── setup.R
│   ├── simulate_from_dag.R
│   ├── test.R
│   └── wald_test.R
├── README.md
├── runEx.sh                              # script for running example.R
├── runSim.sh                             # script for running simulation.R
├── run_simulation.R                      # code for the simulation study
├── sc_graph.R                            # code for the real data example
└── simulation.R                          # wrapper of run_simulation.R
```
