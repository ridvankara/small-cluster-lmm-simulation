# Small-Cluster Inference in Linear Mixed Models

Monte Carlo simulation study comparing REML, Kenward-Roger, 
Satterthwaite, and Parametric Bootstrap.

## Paper
Kara, R. (under review). Small-Cluster Inference in Linear Mixed Models...
*Communications in Statistics — Simulation and Computation*

## Scripts
- `run_simulation.R` — Primary frequentist simulation (360 conditions × 1,000 rep)
- `run_bootstrap.R` — Parametric bootstrap (J=3,5, 72 conditions × 1,000 rep)
- `run_bayes.R` — Preliminary Bayesian analysis (J=3, n=50, 200 rep)

## Requirements
R 4.5.2, lme4 1.1.38, lmerTest 3.2.1, pbkrtest 0.5.5, brms 2.23.0, 
furrr 0.3.1, seed = 2025
