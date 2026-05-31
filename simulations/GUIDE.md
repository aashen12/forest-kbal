# Simulations Guide

This folder runs Monte Carlo simulations for the forest kernel balancing weights paper (Shen et al. 2025). There are four simulation studies, each with its own R script and SLURM submission script.

All simulation scripts share helper functions from the `functions/` folder one level up. You do **not** run those function files directly -- they are sourced automatically.

---

## Quick Start

All commands below assume you are in the `simulations/` directory. On the SLURM cluster:

```bash
cd simulations/

# Simulation 1 (main paper, Section 4): 1000 reps
sbatch sim1.sh 1000

# Simulation 2 (appendix, moderate overlap): 1000 reps
sbatch sim2.sh 30 1000

# Simulation 2 (appendix, strong overlap): 1000 reps
sbatch sim2.sh 100 1000

# Varsweep (sensitivity analysis for mixing parameter p): 2000 reps
sbatch varsweep.sh 2000

# Simulation 1 revision (dimensionality sweep): one job per q, 1000 reps each
for q in 10 40 100; do sbatch sim1-revision.sh $q 1000; done
```

For local testing with fewer reps (no SLURM needed):

```bash
Rscript simulation1.R 5
Rscript simulation2.R 30 5
Rscript varsweep.R 5
Rscript sim1-revision.R 40 5
```

When no SLURM environment is detected, the scripts default to 4 cores.

Once results exist, turn them into the paper figures (no SLURM needed -- see "Generating Figures" below):

```bash
Rscript visualize-sim1.R           # Figure 3 + B.1-B.4 from one sim1 run
Rscript visualize-sim1-revision.R  # dimensionality-sweep figures (overlays q)
```

---

## Simulation 1: Tarr & Imai (2025) DGP

**Script:** `simulation1.R`
**Shell wrapper:** `sim1.sh` (16 cores, `high` partition)

This is the main simulation in Section 4. Generates data with 10 nonlinearly transformed covariates, logistic treatment assignment, and heterogeneous treatment effects. True ATT varies by draw.

**Usage:**

```bash
sbatch sim1.sh <num_reps>      # on SLURM
Rscript simulation1.R <num_reps>   # locally
```

**Parameters:** n = 1000 observations per rep. Pilot sample is half the controls (~250 units).

**Output:** `results/simulation1-<num_reps>-<date>.RData`

**Logs:** `logs/sim1-<date>.txt`

---

## Simulation 2: Kim et al. (2024) DGP

**Script:** `simulation2.R`
**Shell wrapper:** `sim2.sh` (11 cores, `high` partition)

Appendix simulation with 6 covariates (5 normal + 1 binary) and adjustable overlap between treatment and control groups. True ATT = 0 by construction.

**Usage:**

```bash
sbatch sim2.sh <overlap> <num_reps>
Rscript simulation2.R <overlap> <num_reps>
```

The `overlap` parameter controls treatment-control separation:
- `30` = moderate overlap (harder estimation problem)
- `100` = strong overlap (easier estimation problem)

The shell script validates that overlap is either 30 or 100.

**Parameters:** n = 1000 per rep.

**Output:** `results/simulation2-overlap-<overlap>-numsim-<num_reps>-<date>.RData`

**Logs:** `logs/sim2-<overlap>-<date>.txt`

---

## Simulation 1 Revision: Covariate Dimensionality

**Script:** `sim1-revision.R`
**Shell wrapper:** `sim1-revision.sh` (16 cores, `high` partition)

Added in revision to answer Reviewer 1, Comment 3: *how does the dimensionality q of the covariates X affect performance, and how large q can the method handle?* This re-runs Simulation 1 while varying the number of covariates q.

The DGP (`make_data_sim1_highdim` in `../functions/dgp.R`) generalizes the Simulation 1 DGP: it builds q latent Gaussians in blocks of 10 and cycles both the nonlinear covariate transformations and the treatment/outcome model over each block, so the added covariates are genuine confounders (4 active + 6 noise per block). The propensity and outcome linear predictors are rescaled by sqrt(number of blocks), which holds propensity overlap and outcome signal-to-noise constant as q grows. This isolates the effect of dimensionality. At q = 10 the DGP reduces exactly to Simulation 1.

**Usage:**

```bash
sbatch sim1-revision.sh <q> <num_reps>
Rscript sim1-revision.R <q> <num_reps>
```

The `q` parameter is the number of covariates; the shell script validates it is one of `{10, 40, 100}` (all multiples of 10). To sweep all dimensions, submit one job per q:

```bash
for q in 10 40 100; do sbatch sim1-revision.sh $q 1000; done
```

To add more dimensions, append values to the `PARAMS` array in `sim1-revision.sh` (the DGP itself accepts any q, including non-multiples of 10, in which case the final block is partial).

**Parameters:** n = 1000 per rep (fixed across q, so the q/n ratio grows with q).

**Output:** `results/sim1-revision-q-<q>-numsim-<num_reps>-<date>.RData`. Each results row is tagged with a `q` column so the per-q files can be row-bound and compared directly.

**Logs:** `logs/sim1-revision-q<q>-<date>.txt`

---

## Varsweep: Mixing Parameter Sensitivity

**Script:** `varsweep.R`
**Shell wrapper:** `varsweep.sh` (21 cores, `epurdom` partition)

Sweeps the mixing parameter p from 0 to 1 in increments of 0.05, where p controls the weight between raw covariates and kernel features: X_combined = (1-p) * X_raw + p * X_kernel. Uses the same DGP as Simulation 1.

**Usage:**

```bash
sbatch varsweep.sh <num_reps>
Rscript varsweep.R <num_reps>
```

**Parameters:** n = 1000 per rep, p in {0.00, 0.05, 0.10, ..., 1.00}.

**Output:** `results/varsweep-<num_reps>-<date>.RData`

**Logs:** `logs/varsweep-<date>.txt`

---

## What Each Estimator Produces

Every simulation rep evaluates these methods across multiple feature representations (raw, rf_only, rf_plus, bart_only, bart_plus, kbal_only, kbal_plus, and full-kernel variants):

| Estimator | Column name |
|-----------|------------|
| L2 Balancing Weights | `bal.wgt` |
| Inverse Probability Weights | `ipw` |
| OLS Outcome Regression | `ols` |
| Augmented L2 BW (doubly robust) | `ols_l2` |
| Augmented IPW (doubly robust) | `ols_ipw` |

For each estimator, the output includes: estimated ATT, SE, 95% CI, bias, coverage, percent bias reduction (PBR), and effective sample size (ESS).

---

## Folder Structure

```
simulations/
  simulation1.R        # Main simulation script (Section 4)
  simulation2.R        # Appendix simulation script
  sim1-revision.R      # Dimensionality sweep (revision, reviewer comment 3)
  varsweep.R           # Mixing parameter sensitivity analysis
  sim1.sh              # SLURM wrapper for simulation1.R
  sim2.sh              # SLURM wrapper for simulation2.R
  sim1-revision.sh     # SLURM wrapper for sim1-revision.R
  varsweep.sh          # SLURM wrapper for varsweep.R
  visualize-sim1.R          # Reproduce Figure 3 + B.1-B.4 from a sim1 run
  visualize-sim1-revision.R # Dimensionality-sweep figures (overlays q)
  results/             # .RData output files
  logs/                # Timestamped progress logs
  paper-figs/          # PDF figures for the paper
  slurm_output/        # SLURM stdout/stderr

functions/  (one level up, sourced automatically)
  dgp.R                # DGP definitions: make_data_sim1(), make_data_sim1_highdim(), make_data_sim2()
  sim-utils.R          # Shared utilities: logging, parallel setup, result processing
  sim-eval-funcs.R     # Main evaluation pipeline: eval_data()
  sim-estimation-funcs.R  # Core estimators: balancingWeights(), logisticIPW(), etc.
  BART-features.R      # BART kernel extraction and PCA
  randomForestFeatures.R  # Random Forest kernel extraction and PCA
  varsweep-funcs.R     # Varsweep-specific evaluation: eval_data_varsweep()
  sim-plot-utils.R     # Figure helpers: summarise_sim_results(), plot_sim_metric()
```

---

## Shared Function Files (do not run directly)

These live in `../functions/` and are sourced by the simulation scripts:

| File | What it provides |
|------|-----------------|
| `dgp.R` | `make_data_sim1(n)`, `make_data_sim1_highdim(n, q)`, and `make_data_sim2(n, sig.ep)` -- generate simulated datasets |
| `sim-utils.R` | Package loading, logging, parallel setup, result flattening |
| `sim-eval-funcs.R` | `eval_data()` -- the main workhorse that generates kernel features and runs all estimators |
| `sim-estimation-funcs.R` | Individual estimator functions (balancing weights, IPW, OLS, augmented) |
| `BART-features.R` | `bart_kernel_matrix()` and `pca_bart()` -- extract BART implied kernel |
| `randomForestFeatures.R` | `rf_kernel_matrix()` -- extract RF implied kernel |
| `varsweep-funcs.R` | `eval_data_varsweep()` -- evaluation with mixing parameter p |
| `sim-plot-utils.R` | `summarise_sim_results()`, `plot_sim_metric()` (single run), and `plot_sim_metric_byq()` (revision: faceted by q) -- shared figure helpers (requires tidyverse) |

---

## Generating Figures

After a simulation has produced results, two scripts turn them into the paper PDFs (written to `paper-figs/`). These are plain `Rscript` jobs -- no SLURM needed -- and only require `tidyverse`.

### `visualize-sim1.R` -- main paper figures

Reproduces the Simulation 1 figures from a single results file. Edit `results_file` at the top to point at the run you want, then:

```bash
Rscript visualize-sim1.R
```

| Figure | Content | Output |
|--------|---------|--------|
| Figure 3 (main) | Relative RMSE, balancing weights | `sim1_relrmse.pdf` |
| Figure B.1 | Absolute relative bias, balancing weights | `sim1_relbias.pdf` |
| Figure B.2 | Percent bias reduction; effective sample size | `sim1_pbr.pdf`, `sim1_ess.pdf` |
| Figure B.3 | Absolute relative bias, trimmed IPW | `sim1_relbias_ipw.pdf` |
| Figure B.4 | Relative RMSE, trimmed IPW | `sim1_relrmse_ipw.pdf` |

Each figure facets Kernel Only vs. Kernel + Raw, with one line per kernel family (Design Kernel = red, RF = blue, BART = orange) and a dashed raw-covariate baseline.

### `visualize-sim1-revision.R` -- dimensionality-sweep figures

Shows how performance degrades as the covariate dimension q grows. Each metric produces **two figures** -- one for **Kernel + Raw** and one for **Kernel Only** -- and each figure is **faceted by q** (one panel per dimension, left to right). Within a panel, color/shape encode the kernel family (Design Kernel, RF, BART) and a thin grey dashed line marks the raw-covariate baseline at that q. The y-axis is shared across panels for direct comparison. (This replaces the earlier single-figure overlay, which packed 9 lines into each panel.) The relative-RMSE pair is the headline: `sim1revision_relrmse_kernelraw.pdf` is the main-text figure and `sim1revision_relrmse_kernelonly.pdf` the appendix version.

List the files to plot in the `results_files` vector at the top of the script. It is a *named* vector: each name is the covariate dimension q and each value is the file path. The loader stamps every file with its q, so the q = 10 publication file (`results-revision/simulation1-1000-Nov-12-2025.RData` -- just the Simulation 1 output, with no q column of its own) is baked in alongside the `sim1-revision-q-*` files. Add or comment out lines as runs finish.

```bash
Rscript visualize-sim1-revision.R
```

Each metric is written as `sim1revision_<metric>_kernelraw.pdf` and `sim1revision_<metric>_kernelonly.pdf`. If a panel looks compressed because a few-PC point spikes, pass an explicit `ylim` to `plot_sim_metric_byq()` (it zooms via `coord_cartesian`, so no data is dropped).

---

## Notes

- Simulation 1 results are used in the main paper (Section 4). Simulation 2 and varsweep are in the appendix. The Simulation 1 revision (dimensionality sweep) was added to address Reviewer 1, Comment 3.
- Each simulation writes progress to a log file in `logs/`. Check these if a run seems stuck.
- The number of principal components tested varies: Simulation 1 (and its revision) uses nc in {2, 5, 10, 15, 25, 50, 100}; Simulation 2 and WSC/soldiering analyses use nc in 2:50. Note that nc (kernel principal components) is distinct from q (covariate dimension) swept in the revision.
- To reproduce paper results, use the same number of reps and let the default seeds run. Seeds are deterministic given the rep count.
- The revision's `sim1-revision.R` at q = 10 reproduces Simulation 1 exactly (same DGP), so it can double as a consistency check.
