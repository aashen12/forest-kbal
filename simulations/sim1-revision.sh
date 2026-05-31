#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%j.o
#SBATCH --error=slurm_output/%x_%j.e
#SBATCH -c 30
#SBATCH -p epurdom

# Simulation 1 revision: covariate dimensionality experiment (reviewer comment 3)
# Usage: sbatch sim1-revision.sh <q> <num_reps>
# Example: sbatch sim1-revision.sh 40 1000
#
# Submit one job per dimension to sweep q in parallel on the cluster:
#   for q in 10 40 100; do sbatch sim1-revision.sh $q 1000; done

# Work from this script's directory. Under SLURM, $0 is the spooled copy of the
# script, so prefer SLURM_SUBMIT_DIR (the directory you ran sbatch from).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Valid covariate dimensions (multiples of 10). To add more, append values here.
PARAMS=(10 40 50 100)

q=$1
num_sim=$2

# Validate q parameter
VALID=false
for v in "${PARAMS[@]}"; do
  if [[ "$q" == "$v" ]]; then
    VALID=true
    break
  fi
done

if ! $VALID; then
  echo "ERROR: invalid q '$q' (must be one of: ${PARAMS[*]})"
  exit 1
fi

echo "Running sim1-revision.R with q=$q, num_sim=$num_sim"
Rscript sim1-revision.R "$q" "$num_sim"
