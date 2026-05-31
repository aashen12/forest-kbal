#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aashen@berkeley.edu
#SBATCH --output=slurm_output/%x_%A_%a.o
#SBATCH --error=slurm_output/%x_%A_%a.e
#SBATCH -c 2
#SBATCH --mem-per-cpu=24G
#SBATCH -p high
#SBATCH --array=1-100

# Simulation 1 revision via a SLURM job array. Each array task runs
# <reps_per_job> reps of the dimensionality-sweep DGP with a unique seed and
# writes results-spawn/sim1-revision-q-<q>-n-<n>-task-<id>.RData. After all
# tasks finish, combine them with combine-spawn-results.R.
#
# Total reps = (array size) * reps_per_job.  Default below = 100 * 10 = 1000.
#
# Args:
#   q             Number of covariates (e.g. 10, 50, 100)
#   reps_per_job  Monte Carlo reps run by each array task
#   n             Optional sample size (default 1000). Pass 5000 for the n=5000
#                 sweep requested by the collaborators.
#
# Usage:
#   sbatch sim1-revision-array.sh <q> <reps_per_job> [n]
# Examples:
#   sbatch sim1-revision-array.sh 50 10                       # 100 tasks x 10 reps = 1000 reps at n=1000, q=50
#   sbatch sim1-revision-array.sh 100 10                      # same, q=100
#   sbatch sim1-revision-array.sh 50 10 5000                  # SAME structure, but n=5000 (slower per rep)
#   sbatch sim1-revision-array.sh 100 10 5000                 # n=5000, q=100
#   sbatch --array=1-200 sim1-revision-array.sh 50 5 5000     # finer split: 200 x 5 reps at n=5000
#   sbatch --array=1-100%20 sim1-revision-array.sh 50 10 5000 # cap to 20 concurrent tasks

# Under SLURM, $0 is the spooled copy; prefer SLURM_SUBMIT_DIR.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Keep BLAS/OMP single-threaded so many array tasks on one node don't
# oversubscribe its cores. The within-task parallelism comes from mclapply,
# which respects SLURM_CPUS_PER_TASK (-c above).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

q=$1
reps_per_job=$2
n=${3:-1000}     # optional sample size, default 1000

if [[ -z "$q" || -z "$reps_per_job" ]]; then
  echo "ERROR: usage: sbatch sim1-revision-array.sh <q> <reps_per_job> [n]"
  exit 1
fi

echo "Array task ${SLURM_ARRAY_TASK_ID:-(local)}: q=${q}, reps_per_job=${reps_per_job}, n=${n}"
Rscript sim1-revision-array.R "$q" "$reps_per_job" "$n"
