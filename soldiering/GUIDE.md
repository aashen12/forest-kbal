# Soldiering Application Guide

This folder contains the empirical analysis of the Blattman & Annan (2010) child soldiering study in Uganda (Section 5.2 of Shen et al. 2025). The outcome is years of education. Treatment is abduction by the Lord's Resistance Army (`abd`).

---

## Quick Start

### Running the main analysis (SLURM cluster)

```bash
cd soldiering/
sbatch soldiering-xfit.sh
```

### Running locally (slow -- takes many hours with 10 repeats)

```bash
cd soldiering/
Rscript soldiering-xfit.R
```

### Generating paper figures (after results exist)

```bash
cd soldiering/
Rscript visualize.R
```

---

## Scripts: What to Run

### 1. `soldiering-xfit.R` -- Main Cross-Fit Analysis

This is the primary analysis script. It implements Algorithm 1 (2-fold cross-fitting) with 10 random splits of the control group.

**What it does:**
1. Loads `blattman_claude.csv` and defines covariates via `soldiering-setup.R`
2. Splits controls 50/50 into two halves
3. Each half takes a turn as the pilot (training) sample
4. Estimates ATT using all methods (balancing weights, IPW, OLS, augmented)
5. Averages results across folds and repeats

**How to run:**

```bash
sbatch soldiering-xfit.sh    # on SLURM (11 cores, yss partition)
Rscript soldiering-xfit.R    # locally (defaults to 4 cores)
```

**Output:** `results/soldiering-xfit-educ-2010.RData`

**Logs:** `logs/soldiering-xfit-<date>.txt`

**Key parameters:**
- `n_repeat = 10` (10 random control splits)
- `seeder = 2011` (determines output filename and seed sequence)
- Outcome: years of education (`educ`)
- Parallel-safe RNG: L'Ecuyer-CMRG

---

### 2. `pilot-data-analysis.R` -- Double-Dip Baseline (No Cross-Fitting)

Uses ALL controls as both the pilot (training) sample and part of the analysis sample. This is the "double-dip" baseline where the same data appears in both roles -- no sample splitting.

**What it does:**
1. Loads data and removes covariates with missingness in the control group
2. Runs a single `eval_data()` call with `dataset = "soldiering-dbldip"`

**How to run:**

```bash
Rscript pilot-data-analysis.R
```

**Output:** `results/multi-dataset-results-educ-dbldip.RData`

---

### 3. `visualize.R` -- Paper Figures (Cross-Fit Results)

Post-processes results from `soldiering-xfit.R` and creates the figures for the paper.

**How to run:**

```bash
Rscript visualize.R
```

**Input:** `results/soldiering-xfit-educ-2010.RData`

**Outputs:**
| File | Description |
|------|-------------|
| `paper-figs/soldiering_full_plot.pdf` | ATT estimates by number of principal components (main result) |
| `paper-figs/soldiering_tf_plot.pdf` | Forest plot at nc = 5 (point estimates + 95% CIs) |
| `paper-figs/soldiering_scree_plot.pdf` | Scree plots for RF and BART kernels across repeats |

Also computes (but does not save) ESS/PBR metrics plots.

---

### 4. `pilot-data-vis.R` -- Pilot Data Figures (Double-Dip Results)

Visualizes results from `pilot-data-analysis.R`. Similar plots to `visualize.R` but for the no-cross-fitting baseline.

**How to run:**

```bash
Rscript pilot-data-vis.R
```

**Input:** `results/multi-dataset-results-educ-dbldip.RData`

---

### 5. `create_unfound_ctrl.R` -- Data Preparation (One-Time)

Extracts unfound controls from the full Blattman dataset. Only needs to be run once to create `unfound_ctrl.csv`.

**How to run:**

```bash
Rscript create_unfound_ctrl.R
```

**Input:** `consequences_analysis_full.csv`

**Output:** `unfound_ctrl.csv`

---

## Scripts: Sourced by Other Scripts (Do Not Run Directly)

### `soldiering-setup.R`

Defines two functions used by both analysis scripts:
- `load_soldiering_data(path)` -- reads `blattman_claude.csv`, renames `abd` to `Z` and `educ` to `Y`
- `define_soldiering_covariates(raw.data)` -- returns named list of covariate groups:
  - `AC_vars`: age and geographic region indicators
  - `ACG_vars`: AC_vars + war experience indicators (G1-G8)
  - `ctrls`: household and parental controls with polynomial expansions
  - `ctrl_hh`: household-level controls (raw)
  - `ctrl_i`: individual parental/orphan controls
  - `full_covs`: union of all groups

### Shared functions in `../functions/`

| File | Key functions |
|------|--------------|
| `cross-fit.R` | `run_cross_fit()`, `process_applied_results()` |
| `sim-utils.R` | `load_sim_packages()`, `setup_parallel()`, `create_log()` |
| `sim-eval-funcs.R` | `eval_data()` -- main evaluation pipeline |
| `sim-estimation-funcs.R` | `balancingWeights()`, `logisticIPW()`, etc. |
| `BART-features.R` | `bart_kernel_matrix()`, `pca_bart()` |
| `randomForestFeatures.R` | `rf_kernel_matrix()` |

---

## Folder Structure

```
soldiering/
  soldiering-xfit.R       # Main analysis (RUN THIS)
  soldiering-xfit.sh      # SLURM wrapper for above
  soldiering-setup.R      # Data loading + covariate definitions (sourced)
  pilot-data-analysis.R   # Double-dip baseline analysis
  visualize.R             # Paper figure generation
  pilot-data-vis.R        # Double-dip figure generation
  create_unfound_ctrl.R   # One-time data prep
  blattman_claude.csv     # Main dataset
  consequences_analysis_full.csv  # Full Blattman dataset
  unfound_ctrl.csv        # Unfound controls (created by create_unfound_ctrl.R)
  results/                # .RData output files
  logs/                   # Timestamped progress logs
  paper-figs/             # PDF figures for the paper
  slurm_output/           # SLURM stdout/stderr
```

---

## Typical Workflow

1. Make sure `blattman_claude.csv` is in this folder
2. Submit the cross-fit analysis: `sbatch soldiering-xfit.sh`
3. Wait for it to finish (check `logs/` for progress, `slurm_output/` for errors)
4. Generate figures: `Rscript visualize.R`
5. Figures appear in `paper-figs/`

For the double-dip baseline, run `pilot-data-analysis.R` then `pilot-data-vis.R`.

---

## Notes

- The cross-fit analysis takes several hours on 10 cores. Each of the 10 repeats runs two full `eval_data()` calls (one per fold), and each call generates RF, BART, and KBal kernels.
- The `seeder` value (2011) determines the output filename (`soldiering-xfit-educ-2010.RData`). If you change it, update `visualize.R` to load the correct file.
- Results in `results/` include several older runs with different seeds and outcomes (distress, logged education). The paper uses `soldiering-xfit-educ-2010.RData`.
