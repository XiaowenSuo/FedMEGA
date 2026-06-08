# Simulation scripts

This directory contains the manuscript simulation scripts for the three-client
FedMEGA analyses.

- `run_binary_simulation_15w20w.R`: binary-trait FedBMM simulation script for
  the 15000samples/20wSNPs setting.
- `run_continuous_simulation_15w20w.R`: continuous-trait FedLMM simulation
  script for the 15000samples/20wSNPs setting.
- `three_site_usage_template.R`: compact template showing how the cleaned core
  functions can be called with three site-specific input files.

The scripts do not include protected genotype or phenotype data. File names and
covariate names should be edited to match the local analysis setting.
