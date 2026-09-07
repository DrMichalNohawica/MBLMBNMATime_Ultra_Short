# Ultra-short vs longer dental implants: marginal bone loss NMA

Analysis code and data for the Bayesian network meta-analysis of marginal
bone loss around ultra-short versus longer dental implants
(manuscript DMP-04508, *Dental and Medical Problems*).

PROSPERO registration: CRD42023400933 (PDF included).

## Contents

- `data/MBL_data_extracted.csv` — the model input: one row per
  study × treatment group × timepoint (months from prosthetic loading),
  with implants (n), mean marginal bone loss (mm) and SD.
- `data/extended_regression.csv` — arm-level table with implant geometry
  (nominal and endosseous length, diameter, platform level, machined
  collar height) and measurement metadata.
- `data/implant_systems.csv` — the per-system catalogue of length
  conventions (whether stated length includes the transmucosal collar),
  verified against manufacturer engineering documents; sources cited
  per row.
- `data/rob_assessments.csv` — risk-of-bias judgements (RoB 2 /
  ROBINS-I) per modelled study.
- `model/multilevel_model.stan` — the longitudinal meta-regression model.
- `scripts/` — the analyses, numbered in run order.

## Run order

```
Rscript scripts/01_netmeta.R            # frequentist per-timepoint NMA + funnels
Rscript scripts/02_mbnmatime.R          # MBNMAtime time-course model (primary)
Rscript scripts/03_stan_regression.R    # Stan longitudinal meta-regression
Rscript scripts/04_contrasts.R          # contrasts + equivalence probabilities
Rscript scripts/06_sensitivity_4mm.R    # <=4 mm (registered criterion) analysis
Rscript scripts/07_bayes_4mm.R          # Bayesian one-timepoint <=4 mm model
Rscript scripts/11_league_tables.R      # full pairwise league tables
```

Outputs are written to `output/`; fitted model objects to `fits/`
(both git-ignored). Scripts 02 and 03 sample for some tens of minutes.

## R dependencies

dplyr, readr, tidyr, stringr, meta, netmeta, MBNMAtime, rstan, metafor.

## Licence / citation

Data tables are provided for reproducibility and audit of the published
analysis; please cite the manuscript.
