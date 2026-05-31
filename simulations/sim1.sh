#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 16
#SBATCH -p high

# Simulation 1: Tarr & Imai (2025) DGP — main paper simulation
# Usage: sbatch sim1.sh <num_reps>
# Example: sbatch sim1.sh 1000

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript simulation1.R $1
