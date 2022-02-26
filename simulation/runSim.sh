#!/bin/bash
#SBATCH --job-name=runSim
#SBATCH --array=1-200
#SBATCH --output=output/runSim-%a.out

Rscript slow-sim.R
