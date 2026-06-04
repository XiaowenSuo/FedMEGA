library(Rcpp)

cpp_files <- setdiff(
  list.files(file.path("code", "cpp"), pattern = "\\.cpp$", full.names = TRUE),
  file.path("code", "cpp", "RcppExports.cpp")
)

for (cpp in cpp_files) {
  sourceCpp(cpp)
}

r_files <- setdiff(
  list.files(file.path("code", "R"), pattern = "\\.R$", full.names = TRUE),
  file.path("code", "R", "RcppExports.R")
)

for (r_file in r_files) {
  source(r_file)
}

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

covariate_names <- c("sex", "age", "PC1", "PC2", "PC3")

fit_binary <- fed_bmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = covariate_names,
  snp_file = "snp200000.txt",
  n_ratio_snp = 30
)

head(fit_binary$assoc)

fit_continuous <- fed_lmm_three_site(
  pheno_files = pheno_files,
  bed_files = bed_files,
  phenotype = "pheno",
  covariates = covariate_names,
  snp_file = "snp200000.txt",
  n_ratio_snp = 100
)

head(fit_continuous$assoc)
