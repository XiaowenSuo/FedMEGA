# FedMEGA

FedMEGA is the core R/C++ implementation of the federated mixed-model GWAS
framework described in the manuscript. The current release implements the
three-client federated setting used in the study.

## Scope

This repository provides the algorithmic implementation for:

- FedBMM: binary-trait federated generalized linear mixed-model GWAS.
- FedLMM: continuous-trait federated linear mixed-model GWAS.
- Federated null regression fitting.
- Federated null mixed-model fitting.
- Federated score testing with GRAMMAR-Gamma approximation.

The current implementation is fixed to three logical clients named `N`, `M`,
and `S`, matching the manuscript analyses. Extending the implementation to an
arbitrary number of sites requires adapting the client aggregation steps.

This repository does not include protected genotype, phenotype, or UK Biobank
individual-level data. Users must supply their own site-specific phenotype,
covariate, and genotype files.

## Installation

```r
install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USER/FedMEGA")
```

For local installation during review:

```r
remotes::install_local("path/to/FedMEGA")
```

## Input contract

Phenotype files should contain:

- one phenotype column, default `pheno`;
- covariate columns specified by the user.

The manuscript simulation examples used `Gender` and the first five genetic
principal components. This is not an algorithmic restriction; any set of
covariates can be supplied through the `covariates` argument as long as the
same columns are available at all three sites.

BED files should be supplied as PLINK `.bed` prefixes accepted by
`BEDMatrix::BEDMatrix()`.

All vectors must be named with the three site identifiers:

```r
pheno_files <- c(
  N = "inputN15w20w.txt",
  M = "inputM15w20w.txt",
  S = "inputS15w20w.txt"
)

bed_files <- c(
  N = "sampleN15w20w.bed",
  M = "sampleM15w20w.bed",
  S = "sampleS15w20w.bed"
)
```

## Limitations of this release

- The public core implementation is restricted to three sites.
- No protected genotype or phenotype data are distributed with the package.
- The binary-trait GRAMMAR-Gamma ratio implementation currently uses an explicit
  projection matrix for the ratio-estimation step. This is appropriate for the
  manuscript-scale simulation setting, but large empirical-scale deployments
  should use a memory-optimized implementation that avoids storing the full
  projection matrix.
- The package is intended for GitHub/Zenodo distribution as manuscript code. It
  is not yet submitted to CRAN.

## Code availability statement

This repository provides the core three-client FedMEGA implementation used in
the manuscript, including FedLMM, FedBMM, federated null-model fitting,
federated null mixed-model fitting, score testing, and GRAMMAR-Gamma
approximation. Protected UK Biobank data and scripts requiring access to
individual-level UK Biobank records are not included. Users must obtain access
to protected datasets through the corresponding data access procedures.

## Basic usage

Binary trait:

```r
library(FedMEGA)

fit_bmm <- fed_bmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = c("Gender", "PC1", "PC2", "PC3", "PC4", "PC5"),
  snp_file = "snp200000.txt",
  n_ratio_snp = 30
)

head(fit_bmm$assoc)
```

Continuous trait:

```r
fit_lmm <- fed_lmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = c("Gender", "PC1", "PC2", "PC3", "PC4", "PC5"),
  snp_file = "snp200000.txt",
  n_ratio_snp = 100
)

head(fit_lmm$assoc)
```

## Manuscript scripts

The package exposes cleaned core functions rather than the original
manuscript-run scripts. The original scripts contained local paths and
protected-data assumptions and are therefore not distributed as package APIs.
