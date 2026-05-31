#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 11
#SBATCH -p yss

# Soldiering cross-fit analysis: Blattman & Annan (2010), outcome = education
# Usage: sbatch soldiering-xfit.sh

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript soldiering-xfit.R
