# globalgender

This repository contains data and analysis code to reproduce results from the article **"Demographic Processes Constrain Global Growth in Gender Egalitarianism"**.

### Description

All processed data are stored in `/data`. To reproduce data processing (conducted in `/R/00_data.R`), original data files from external sources (see `/data/raw/readme.md`) need to be placed inside `/data/raw/`.

The scripts stored in `/R` contain the following:

- `functions.R`: Custom helper functions
- `00_data.R`: Cleaning, harmonization and integration of external data sources; export of processed data to `/data`
- `01_hgam.R`, `02_lmm`, `03_asfr.R`, `04_measure.R`: Statistical analysis; export of results to `/results`

The script `/R/03_asfr.R` further imports the Stan script that is stored in `/Stan`.

The Supplementary Materials of the article can be reproduced via the RMarkdown script in the `/report` folder. This requires that some intermediate results (model objects from time-consuming fitting routines) are placed inside `/results/intermediate` by running the analysis scripts in `/R`. The large files are not part of the repository, but they can be recreated.

### Version info

R session info can be accessed <a href="sessionInfo.txt">here</a>. MCMC sampling was conducted via the `cmdstan` backend `v2.37.0`.
