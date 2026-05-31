#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 11
#SBATCH -p yss

# WSC cross-fit analysis: Byrd (2021) within-study comparison
# Usage: sbatch wsc-xfit.sh

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript wsc-xfit.R
