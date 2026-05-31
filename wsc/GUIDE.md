# WSC Application Guide

This folder contains the within-study comparison (WSC) empirical analysis using Byrd (2021) data (Section 5.1 of Shen et al. 2025). There are two outcomes: **math** test scores and **vocabulary** test scores. Each outcome has an experimental benchmark ATT (math = 0.79, vocab = 2.18) that the observational estimators should recover.

---

## Quick Start

### Running the main analysis (SLURM cluster)

```bash
cd wsc/
sbatch wsc-xfit.sh
```

### Running locally

```bash
cd wsc/
Rscript wsc-xfit.R
```

### Generating paper figures (after results exist)

```bash
cd wsc/
Rscript visualize-xfit-paper.R
```

---

## Scripts: What to Run

### 1. `wsc-xfit.R` -- Main Cross-Fit Analysis

This is the primary analysis script. It implements Algorithm 1 (2-fold cross-fitting) with 10 random splits of the control group. Unlike the soldiering analysis, each repeat runs **two** cross-fits: one with untransformed covariates and one with log-transformed covariates.

**What it does:**
1. Loads WSC observational data and defines 25 covariates via `wsc-setup.R`
2. Creates a log-transformed copy of the data: log(x + 1) for continuous non-negative covariates
3. Splits controls 50/50 into two halves
4. Runs cross-fitting on both untransformed and log-transformed versions
5. Saves results for downstream visualization

**How to run:**

```bash
sbatch wsc-xfit.sh          # on SLURM (11 cores, yss partition)
Rscript wsc-xfit.R          # locally (defaults to 4 cores)
```

**Key parameters to set before running** (edit lines 29-30 of `wsc-xfit.R`):

| Parameter | Options | Default | Effect |
|-----------|---------|---------|--------|
| `outcome` | `"math"` or `"vocab"` | `"math"` | Which test score to analyze |
| `log_trans` | `TRUE` or `FALSE` | `TRUE` | Whether to include log-transformed cross-fit |

**Output files** (depend on parameter settings):

| outcome | log_trans | Output file |
|---------|-----------|-------------|
| math | TRUE | `results/wsc-math-xfit.RData` |
| math | FALSE | `results/wsc-math-xfit-exp.RData` |
| vocab | TRUE | `results/wsc-vocab-xfit.RData` |
| vocab | FALSE | `results/wsc-vocab-xfit-exp.RData` |

**Logs:** `logs/wsc-xfit-<date>.txt`

**To produce all four result files,** you need to run the script four times, changing `outcome` and `log_trans` each time. The paper uses the math outcome with `log_trans = TRUE`.

---

### 2. `visualize-xfit-paper.R` -- Main Paper Figures

The primary visualization script. Produces all WSC figures for the paper.

**How to run:**

```bash
Rscript visualize-xfit-paper.R
```

**Key parameter:** Set `transform` on line 19 to `"log"` or `"exp"` to choose which results file to load.

**Input:** `results/wsc-math-xfit.RData` (log) or `results/wsc-math-xfit-exp.RData` (exp)

**Outputs:**

| File | Description |
|------|-------------|
| `paper-figs/wsc_full_log.pdf` | ATT by number of principal components, faceted by transform |
| `paper-figs/wsc_main.pdf` | Forest plot: all methods, both transforms, at nc = 5 |
| `paper-figs/wsc_main_transformed.pdf` | Forest plot: transformed covariates only |
| `paper-figs/wsc_scree_plot.pdf` | Scree plots for RF and BART kernels |

Also produces PBR/ESS diagnostic plots (displayed but not saved to file).

---

### 3. `visualize-xfit-bstars.R` -- Supplementary Figures

Creates additional visualizations using the exponential transform results.

**How to run:**

```bash
Rscript visualize-xfit-bstars.R
```

**Input:** `results/wsc-math-xfit-exp.RData`

**Outputs:** Line plot by nc (both transforms) and forest plot at nc = 3. Displayed interactively; not saved to file by default.

---

### 4. `visualize-xfit.R` -- Quick Visualization Helper

Wraps plotting in a `create_plot_wsc()` function for easy exploration of different transform levels.

**How to run:**

```bash
Rscript visualize-xfit.R
```

**Input:** `results/wsc-math-xfit-exp.RData`

---

## Scripts: Sourced by Other Scripts (Do Not Run Directly)

### `wsc-setup.R`

Defines two functions:
- `load_wsc_data(outcome)` -- loads `data/wsc-obs-study-data.RData`, sets up Y and Z based on outcome ("math" or "vocab")
- `log_transform_covariates(data.pre)` -- applies log(x + 1) to continuous non-negative covariates with 3+ distinct values; leaves binary variables and Z/Y unchanged

### `Math-cfit.R`

**Legacy file -- do not use.** This is an older version of the analysis that sources from an `R/` subdirectory that no longer exists. It has been superseded by `wsc-xfit.R`.

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
wsc/
  wsc-xfit.R              # Main analysis (RUN THIS)
  wsc-xfit.sh             # SLURM wrapper for above
  wsc-setup.R             # Data loading + covariate definitions (sourced)
  visualize-xfit-paper.R  # Paper figures (RUN THIS after results exist)
  visualize-xfit-bstars.R # Supplementary figures
  visualize-xfit.R        # Quick visualization helper
  Math-cfit.R             # LEGACY -- do not use
  data/
    wsc-obs-study-data.RData  # Byrd (2021) observational data
    covs.RData                # Covariate definitions (legacy, not used by current scripts)
  results/                # .RData output files
  logs/                   # Timestamped progress logs
  paper-figs/             # PDF figures for the paper
  slurm_output/           # SLURM stdout/stderr
```

---

## Covariates (25 total)

The WSC analysis uses these covariates from the Byrd (2021) observational study:

| Category | Variables |
|----------|-----------|
| Demographics | female, white, black, asian, hisp, married, logAge, income |
| Education | collegeS, collegeM, collegeD, calc |
| Personality (Big 5) | big5O, big5C, big5E, big5A, big5N |
| Psych measures | AMAS, logBDI, MCS, GSES |
| Reading habits | logBooks, mathLike |
| Pre-test scores | vocabPre, mathPre |

---

## Typical Workflow

1. Make sure `data/wsc-obs-study-data.RData` exists
2. Set `outcome` and `log_trans` in `wsc-xfit.R` (default: math, log)
3. Submit: `sbatch wsc-xfit.sh`
4. Wait for completion (check `logs/` for progress)
5. Generate figures: `Rscript visualize-xfit-paper.R`
6. Figures appear in `paper-figs/`

To analyze vocabulary scores, change `outcome <- "vocab"` in `wsc-xfit.R` and rerun.

---

## Notes

- The cross-fit analysis takes several hours on 10 cores. Each repeat runs two full cross-fits (untransformed + log-transformed), and each cross-fit involves two `eval_data()` calls.
- The experimental benchmark is shown as a reference line in the plots: 0.79 for math, 2.18 for vocab. These come from the randomized experiment in Byrd (2021).
- `data/covs.RData` is a legacy file from an older version of the analysis. The current scripts define covariates directly in `wsc-setup.R`.
