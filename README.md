# Small-Cluster Inference in Linear Mixed Models

Reproducibility materials for:

> **Small-Cluster Inference in Linear Mixed Models: A Monte Carlo Study of Degrees-of-Freedom Corrections and Bootstrap Procedures**
> Rıdvan Kara — *Communications in Statistics – Simulation and Computation* (under review)

This repository contains the full simulation code that produces every result, table, and figure in the manuscript. The study compares five inference procedures for testing a fixed effect in the Gaussian random-intercept model with few clusters (J = 3–15): the uncorrected asymptotic Wald test, the Satterthwaite and Kenward–Roger degrees-of-freedom corrections, and the parametric and wild cluster bootstraps.

---

## Requirements

- **R** ≥ 4.5.2
- R packages (versions used in the study):
  - `lme4` 1.1.38
  - `lmerTest` 3.2.1
  - `pbkrtest` 0.5.5
  - `lmeresampler` 0.2.4
  - `future`, `furrr` (parallel execution)
  - `dplyr` (summaries)
  - `ggplot2` (figures)

Install with:

```r
install.packages(c("lme4", "lmerTest", "pbkrtest", "lmeresampler",
                   "future", "furrr", "dplyr", "ggplot2"))
```

---

## Repository contents

```
R/
  01_main_simulation.R              Main frequentist design + parametric bootstrap
  02_wild_bootstrap.R               Wild cluster bootstrap (Webb weights)
  03_extensions_covariates_missing.R  Multiple-covariate and missing-data extensions
  04_convergence_sensitivity.R      Convergence / singular-fit rates + exclusion sensitivity
  05_level2_predictor.R             Cluster-level (Level-2) predictor sub-study
  06_make_figures.R                 Generates Figures 1–5 from the result files
results/                            Output .rds files are written here (see Data below)
README.md
```

---

## Run order and outputs

The scripts are independent except that `06_make_figures.R` reads the result
files produced by the others. Run them in the numbered order.

| Script | Produces | Manuscript element |
|---|---|---|
| `01_main_simulation.R` | `sim_freq_results_v2.rds`, `sim_boot_results_v2.rds` | Sections 3.1–3.4; Figures 1–2; Tables 1, 1b |
| `02_wild_bootstrap.R` | `sim_wild_results_v2.rds` | Section 3.5; Figure 3 |
| `03_extensions_covariates_missing.R` | `sim_phase2_results.rds` | Section 3.6; Figures 4–5 |
| `04_convergence_sensitivity.R` | `sim_conv_results.rds` | Section 3.6 (convergence, singular fits, sensitivity) |
| `05_level2_predictor.R` | `sim_level2_results.rds` | Section 3.8; Table 3 |
| `06_make_figures.R` | `figure1.png … figure5.png` | Figures 1–5 |

Each script prints its key summary tables to the console on completion.

---

## Design summary

- **Primary design (degrees-of-freedom procedures):** fully factorial —
  J ∈ {3, 5, 7, 10, 15} × n ∈ {20, 50, 150} × ICC ∈ {0.05, 0.15, 0.30, 0.50}
  × {balanced, unbalanced} × β₁ ∈ {0, 0.3, 0.5} = **360 conditions**, 1,000
  replicates each.
- **Bootstrap procedures (parametric, wild):** prespecified balanced subset
  (n ∈ {20, 50}, ICC ∈ {0.05, 0.30}) across all five cluster counts, 250
  replicates each, B = 499 resamples.
- **Extensions:** multiple covariates and 20% missing data (MCAR/MAR),
  df-based procedures only.
- **Sub-studies:** convergence/singular-fit + exclusion sensitivity (J = 3, 5);
  cluster-level (Level-2) predictor (J = 3–15).

Unbalanced cluster sizes are drawn from `Uniform[0.5n, 1.5n]` (mean = n).
Every condition is seeded reproducibly as `seed = 2025 + condition index`
(with method-specific offsets), so results are exactly reproducible.

---

## Reproducibility notes

- **Paths.** The scripts write to and read from a working directory defined at
  the top of each file (`DIR <- "C:/simulation_lmm"` on the original Windows
  machine). Change `DIR` to any writable folder before running; on
  Linux/macOS use e.g. `DIR <- "~/simulation_lmm"`.
- **Library path.** The `.libPaths(...)` line at the top of each script is
  specific to the original machine and can be removed or edited.
- **Checkpointing.** Long runs save a `*_checkpoint.rds` after every batch, so
  an interrupted run resumes where it stopped. Delete the checkpoint to start
  fresh.
- **Memory.** The bootstrap scripts use `WORKERS = 2` to stay within 8 GB of
  RAM; increase `WORKERS` on machines with more memory.

---

## Data (result files)

Running the scripts regenerates all result files under `results/`. The largest
file (`sim_freq_results_v2.rds`, 360,000 rows × 3 methods) may exceed GitHub's
100 MB per-file limit; if so, host the result files via Git LFS or an archival
data repository (e.g., Zenodo or OSF) and link them here. All figures and
tables in the manuscript can be regenerated from these files with
`06_make_figures.R` and the summary blocks printed by each script.

---

## Citation

If you use this code, please cite the manuscript (details to be updated on
acceptance).

## License

Released under the MIT License.
