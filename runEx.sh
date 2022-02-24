#!/bin/bash
#SBATCH --job-name=runEx
#SBATCH --output=output/runEx.out
#SBATCH --partition=standard

Rscript example.R
