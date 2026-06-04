#' Run FedLMM for a three-site federated continuous-trait GWAS
#'
#' This is the core three-site implementation used by the manuscript scripts.
#'
#' @param pheno_files Named vector with entries `N`, `M`, and `S`.
#' @param bed_files Named vector with entries `N`, `M`, and `S`.
#' @param phenotype Phenotype column name.
#' @param covariates Optional covariate column names. The number and identity of
#'   covariates are user-defined and are not restricted to the manuscript setup.
#' @param snp_file Optional file containing SNP IDs in test order.
#' @param n_ratio_snp Number of SNPs used to estimate the GRAMMAR-Gamma factor.
#' @param maxiter Maximum null mixed-model iterations.
#' @param tol Convergence tolerance for variance-component updates.
#' @param verbose Whether to print fitting progress.
#'
#' @return A list containing `null_model` and `assoc`.
#' @export
fed_lmm_three_site <- function(pheno_files,
                               bed_files,
                               phenotype = "pheno",
                               covariates = NULL,
                               snp_file = NULL,
                               n_ratio_snp = 100,
                               maxiter = 20,
                               tol = 0.02,
                               verbose = TRUE) {
  dat <- read_three_site_phenotypes(pheno_files, phenotype, covariates)
  XN <- dat$N$X
  XM <- dat$M$X
  XS <- dat$S$X
  yN <- dat$N$y
  yM <- dat$M$y
  yS <- dat$S$y

  beta0 <- fed_null_linear_three_site(XN, yN, XM, yM, XS, yS)
  eta0 <- c(drop(XN %*% beta0), drop(XM %*% beta0), drop(XS %*% beta0))

  result <- suoglmmkin.ai_PCG_Rcpp_Quan_no(
    N_n = length(yN),
    N_m = length(yM),
    N_s = length(yS),
    yN = yN,
    yM = yM,
    yS = yS,
    XN = XN,
    XM = XM,
    XS = XS,
    eta0 = eta0,
    alpha0 = beta0,
    isCovariateOffset = FALSE,
    tau = c(0, 0),
    fixtau = c(0, 0),
    tauInit = c(0, 0),
    offset = NULL,
    family = "poisson",
    maxiter = maxiter,
    tol = tol,
    verbose = verbose
  )

  geno <- load_three_site_bed(bed_files)
  snp <- NULL
  if (!is.null(snp_file)) {
    snp <- data.table::fread(snp_file, header = FALSE, data.table = FALSE)[[1]]
  }
  assoc <- score_continuous_three_site(
    result,
    geno,
    yN,
    yM,
    yS,
    snp = snp,
    n_ratio_snp = n_ratio_snp
  )
  list(null_model = result, assoc = assoc)
}
