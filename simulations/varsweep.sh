#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 21
#SBATCH -p epurdom

# Varsweep: Feature weight sensitivity analysis (mixing parameter p)
# Usage: sbatch varsweep.sh <num_reps>
# Example: sbatch varsweep.sh 2000

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript varsweep.R $1
