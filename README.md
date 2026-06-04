# FedMEGA manuscript code

This repository contains the custom analysis code for the FedMEGA manuscript.
FedMEGA is a federated mixed-model GWAS framework with two trait-specific
modules:

- FedBMM for binary phenotypes.
- FedLMM for continuous phenotypes.

The code is provided as a manuscript code archive for GitHub/Zenodo deposition,
not as a CRAN-style R package. It contains the core R/C++ implementation and
simulation scripts used to run the three-client federated analyses described in
the manuscript.

## Repository contents

```text
code/
  R/       Core R functions for FedBMM, FedLMM, null-model fitting, and score tests.
  cpp/     Rcpp/C++ source files used by the core algorithms.

scripts/
  simulation/
    run_binary_simulation.R
    run_continuous_simulation.R
    three_site_usage_template.R
  figures/
    Figure-generation scripts, when available.

environment/
  dependency_versions.tsv

README.md
LICENSE
```

## Scope

The public code implements the three-client horizontal federated setting used in
the manuscript, with logical clients named `N`, `M`, and `S`. Extending the code
to an arbitrary number of participating sites requires adapting the aggregation
steps in the R scripts.

The repository does not include protected genotype, phenotype, or UK Biobank
individual-level data. Users must supply their own site-specific phenotype,
covariate, genotype, and SNP-list files.

## Inputs expected by the simulation scripts

The simulation scripts assume three local clients. Phenotype files should contain
one phenotype column and user-specified covariate columns. In the manuscript
simulation setting, the covariates were sex and five genetic principal
components; this is an analysis choice rather than an algorithmic restriction.

The example file names used by the 15w/20w scripts are:

```text
Binary phenotype files:
  inputN.txt
  inputM.txt
  inputS.txt

Continuous phenotype files:
  inputN_q.txt
  inputM_q.txt
  inputS_q.txt

Genotype BED files:
  sampleN.bed
  sampleM.bed
  sampleS.bed

SNP list:
  snp.txt
```

These file names can be modified directly in the scripts for other simulation
settings.

## Running the simulation scripts

The scripts are intended to be run from a working directory containing the
required phenotype, genotype, SNP-list, and C++ source files. For example:

```r
source("scripts/simulation/run_binary_simulation_15w20w.R")
source("scripts/simulation/run_continuous_simulation_15w20w.R")
```

Because the original analyses were based on protected UK Biobank-derived
genotypes, this repository does not provide runnable toy genotype data.

## Main dependencies

The analyses used R together with `data.table`, `Rcpp`, `RcppArmadillo`,
`BEDMatrix`, `SAIGE`, `pryr`, `psych`, `missMethods`, and base R statistical
functions. A dependency summary is provided in
`environment/dependency_versions.tsv`.

## Code availability statement

Custom code used to implement FedMEGA, including the core FedLMM and FedBMM
algorithms, federated null-model fitting, federated score testing,
GRAMMAR-Gamma approximation, and simulation scripts, is available in this
repository and will be archived with a permanent Zenodo DOI. Protected
individual-level UK Biobank genotype and phenotype data are not distributed and
can be accessed only through the UK Biobank application process.
