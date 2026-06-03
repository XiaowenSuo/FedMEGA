#' Read three-site phenotype and covariate files
#'
#' @param files Named character vector or list with entries `N`, `M`, and `S`.
#' @param phenotype Name of the phenotype column.
#' @param covariates Character vector of covariate column names. The default
#'   matches the manuscript simulation setup; users can provide any covariate
#'   columns shared by all three sites.
#' @param add_intercept Logical; whether to prepend an intercept column.
#'
#' @return A list with entries `N`, `M`, and `S`, each containing `X` and `y`.
#' @export
read_three_site_phenotypes <- function(files,
                                       phenotype = "pheno",
                                       covariates = c("Gender", "PC1", "PC2", "PC3", "PC4", "PC5"),
                                       add_intercept = TRUE) {
  required_sites <- c("N", "M", "S")
  if (is.null(names(files)) || !all(required_sites %in% names(files))) {
    stop("files must be named with entries N, M, and S.")
  }

  read_one <- function(file) {
    dat <- data.table::fread(file, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    missing_cols <- setdiff(c(phenotype, covariates), names(dat))
    if (length(missing_cols) > 0) {
      stop("Missing columns in ", file, ": ", paste(missing_cols, collapse = ", "))
    }
    X <- dat[, covariates, drop = FALSE]
    if (add_intercept) {
      X <- data.frame(Intercept = 1, X)
    }
    list(X = as.matrix(X), y = as.numeric(dat[[phenotype]]))
  }

  list(
    N = read_one(files[["N"]]),
    M = read_one(files[["M"]]),
    S = read_one(files[["S"]])
  )
}

load_three_site_bed <- function(files) {
  required_sites <- c("N", "M", "S")
  if (is.null(names(files)) || !all(required_sites %in% names(files))) {
    stop("bed_files must be named with entries N, M, and S.")
  }

  clean_colnames <- function(geno) {
    new_colnames <- sub("([0-9]+)_.*", "\\1", colnames(geno))
    colnames(geno) <- new_colnames
    geno[, ]
  }

  list(
    N = clean_colnames(BEDMatrix::BEDMatrix(files[["N"]])),
    M = clean_colnames(BEDMatrix::BEDMatrix(files[["M"]])),
    S = clean_colnames(BEDMatrix::BEDMatrix(files[["S"]]))
  )
}
