library(FedMEGA)

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

fit_binary <- fed_bmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = c("Gender", "PC1", "PC2", "PC3", "PC4", "PC5"),
  snp_file = "snp200000.txt",
  n_ratio_snp = 30
)

head(fit_binary$assoc)

fit_continuous <- fed_lmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = c("Gender", "PC1", "PC2", "PC3", "PC4", "PC5"),
  snp_file = "snp200000.txt",
  n_ratio_snp = 100
)

head(fit_continuous$assoc)
