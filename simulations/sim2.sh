#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 11
#SBATCH -p high

# Simulation 2: Kim et al. (2024) DGP with varying overlap — appendix simulation
# Usage: sbatch sim2.sh <overlap> <num_reps>
# Example: sbatch sim2.sh 30 1000

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

PARAMS=(30 100)

overlap=$1
num_sim=$2

# Validate overlap parameter
VALID=false
for v in "${PARAMS[@]}"; do
  if [[ "$overlap" == "$v" ]]; then
    VALID=true
    break
  fi
done

if ! $VALID; then
  echo "ERROR: invalid overlap '$overlap' (must be one of: ${PARAMS[*]})"
  exit 1
fi

echo "Running simulation2.R with overlap=$overlap, num_sim=$num_sim"
Rscript simulation2.R "$overlap" "$num_sim"
