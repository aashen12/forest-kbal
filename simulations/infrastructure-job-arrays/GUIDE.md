# Infrastructure: Job-Array Sweep for Simulation 1 Revision

A faster, cluster-friendly way to run the Simulation 1 revision (the dimensionality sweep, reviewer comment 3). Instead of one SLURM job that runs 1000 reps with `mclapply`, this folder submits a **SLURM job array** that fans out into many small jobs across the cluster.

This folder lives inside `simulations/`. Nothing in the parent `simulations/` folder is modified. If something here breaks, you can always fall back to the original `../sim1-revision.R`.

---

## Why job arrays (and not GPUs)

The bottleneck is BART (`stochtree`), random forests, and `kbal` / SVD — all CPU-bound. None of them are GPU-accelerated in this codebase, so GPUs won't help without rewriting the core estimators. The reps are **embarrassingly parallel** (each one is independent), so a job array is the right tool: 100 independent jobs of 10 reps each run truly in parallel across the cluster's nodes, instead of contending for cores on a single node.

---

## Files in this folder

| File | What it does |
|------|--------------|
| `sim1-revision-array.R` | Array worker. Runs `reps_per_job` reps with a unique seed per rep, writes one task file. |
| `sim1-revision-array.sh` | SLURM array wrapper (`#SBATCH --array=1-100`, 2 cores per task, `high` partition). |
| `combine-spawn-results.R` | Stitches all per-task files into one combined file per q. |
| `results-spawn/` | Output. Holds the per-task `.RData` files and the combined per-q file. |
| `slurm_output/` | SLURM `--output`/`--error` logs, one pair per array task. |

Function files (`../../functions/*.R`) are **referenced**, not duplicated — same DGP and estimation pipeline as `../sim1-revision.R`.

---

## How seeds are handled

Every rep across the whole array gets a unique, reproducible seed of the form

```
seed = SEED_BASE + global_id          # global_id = (task_id - 1) * reps_per_job + local_id
```

`global_id` is unique by construction, so:

- Task 1 reps 1–10 → seeds 23968 … 23977
- Task 2 reps 1–10 → seeds 23978 … 23987
- … (no two reps anywhere ever share a seed → no duplicate data)

The seed is set **explicitly inside each rep** (not via `mc.set.seed`), so the data is reproducible if you rerun.

---

## Workflow

All commands run from inside this folder on the cluster.

### 1. Submit one array per (dimension, sample-size) combination

```bash
cd simulations/infrastructure-job-arrays/

# --- n = 1000 (default; the existing setup) ---
sbatch sim1-revision-array.sh 50  10           # 100 x 10 reps at q=50,  n=1000
sbatch sim1-revision-array.sh 100 10           # same, q=100

# --- n = 5000 (collaborator request: extra row of facets) ---
sbatch sim1-revision-array.sh 50  10 5000      # 100 x 10 reps at q=50,  n=5000
sbatch sim1-revision-array.sh 100 10 5000      # same, q=100, n=5000
```

You don't need to run an array for q=10 at n=1000 — that's just the existing Simulation 1 result (`results-revision/simulation1-1000-Nov-12-2025.RData`), which you can bake into the plot directly. For q=10 at n=5000 there is no pre-existing run, so submit it like the others:

```bash
sbatch sim1-revision-array.sh 10 10 5000       # q=10 at n=5000
```

**Heads-up on n = 5000.** Each rep is roughly an order of magnitude heavier than at n=1000 (kernel matrices are 5000×5000 = 25M entries vs. 1M; BART training on ~2500-row pilots is slower too). Per-task wall time goes from minutes to a few hours. Start with a small smoke test first — e.g. `sbatch --array=1-2 sim1-revision-array.sh 10 1 5000` — and confirm one task finishes cleanly before submitting the full 100-task sweep. You may also want to bump the SLURM time limit with `--time=...` if your partition has a tight default.

**Tuning knobs** (defaults in the `.sh` header):
- `#SBATCH --array=1-100` — number of tasks. Override at submit time: `sbatch --array=1-200 sim1-revision-array.sh 50 5 5000` for 200×5=1000 reps.
- `#SBATCH -c 2` — cores per task. Each task uses `mclapply` over its reps with `SLURM_CPUS_PER_TASK` workers, so bumping `-c` speeds each task up at the cost of running fewer in parallel.
- `#SBATCH --mem-per-cpu=24G` — RAM per allocated CPU (so 48 GB total with `-c 2`). Generous on purpose: at n=5000 kernel matrices are 25M entries (~200 MB) and BART posteriors run hundreds of MB, and 8G was enough to trigger `R_Calloc could not allocate memory` on the heaviest reps. If you still see that error in a task's `.o` log, bump it further.
- `#SBATCH -p high` — partition. Same as the other sim scripts.
- `--array=1-100%20` (mod limit) — cap concurrent tasks if you want to be a polite citizen on the cluster.

Each task writes `results-spawn/sim1-revision-q-<q>-n-<n>-task-<id>.RData`. Logs go to `slurm_output/sim1-revision-array_<jobid>_<taskid>.{o,e}`.

### 2. Watch progress

```bash
squeue -u $USER                                       # see queued/running tasks
ls results-spawn/ | grep task | wc -l                  # how many task files written so far
tail slurm_output/sim1-revision-array_*_1.o           # peek at task 1's output
```

### 3. Combine task files into one data frame per (q, n)

Once all array tasks have finished:

```bash
Rscript combine-spawn-results.R
```

This produces, for each `(q, n)` combination present, a `scenarios` data frame and prints how many task files were folded in. The output naming **preserves your pre-existing n=1000 combined-file names** so nothing in your plotting code has to change:

```
n = 1000 :  results-spawn/sim1-revision-q-<q>-combined-<num_reps>.RData      (legacy name)
n = 5000 :  results-spawn/sim1-revision-q-<q>-n-5000-combined-<num_reps>.RData
```

Legacy task files from before the `n` argument existed (named `sim1-revision-q-<q>-task-*.RData`, no `-n-`) are still picked up and treated as `n = 1000`, so they fold into the legacy-named combined file alongside any new n=1000 task files. The n=1000 and n=5000 sweeps live side by side with no name collision.

### 4. Plot the results

From the `simulations/` folder, run:

```bash
cd ..                          # back up into simulations/
Rscript visualize-sim1-revision.R
```

`visualize-sim1-revision.R` resolves each `(q, n)` cell automatically: you list only the `q` and `n` values in its `results_files` table and it finds the matching combined file in `infrastructure-job-arrays/results-spawn/` by pattern, **ignoring the rep-count suffix** (if several rep counts coexist, the one with the most reps wins). The q = 10 / n = 1000 cell is the lone explicit path (the original Simulation 1 output). Figures are faceted with one **column per q** and one **row per n**, written to `simulations/paper-figs/revision/`. See the header of that script for details.

---

## Local quick test

To sanity-check the worker without SLURM (single task, 1 rep at q=10):

```bash
SLURM_ARRAY_TASK_ID=1 Rscript sim1-revision-array.R 10 1          # n=1000 (default)
SLURM_ARRAY_TASK_ID=1 Rscript sim1-revision-array.R 10 1 5000     # n=5000 (slower)
```

This writes `results-spawn/sim1-revision-q-10-n-<n>-task-0001.RData`. Then `Rscript combine-spawn-results.R` will produce a combined file for each `(q, n)` pair present. Delete those test files when you're done.

---

## What gets pushed to git

The .gitignore is set so:
- Scripts in this folder **are** tracked.
- `slurm_output/` and per-task `.RData` files in `results-spawn/` are **ignored** (high-volume, ephemeral).
- The `.gitkeep` files keep the empty directories around on a fresh clone (so `sbatch` and the worker have somewhere to write).
- **Combined files** (matching `*combined*.RData`) **are** trackable, so you can `git add` and pull them locally for plotting. Or just `scp` them — both work.
