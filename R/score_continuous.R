score_continuous_three_site <- function(result, geno, yN, yM, yS, snp = NULL,
                                        n_ratio_snp = 100, seed = 2025) {
  N_n <- length(yN)
  N_m <- length(yM)
  N_s <- length(yS)
  n <- N_n + N_m + N_s

  sigma <- result$theta[1]
  V_invN <- diag(1 / sigma, N_n)
  V_invM <- diag(1 / sigma, N_m)
  V_invS <- diag(1 / sigma, N_s)

  y0 <- (mean(yN) * N_n + mean(yM) * N_m + mean(yS) * N_s) / n
  yN0 <- yN - y0
  yM0 <- yM - y0
  yS0 <- yS - y0

  mean_N <- colMeans(geno$N)
  mean_M <- colMeans(geno$M)
  mean_S <- colMeans(geno$S)
  global_mean <- (mean_N * N_n + mean_M * N_m + mean_S * N_s) / n

  genoN0 <- sweep(geno$N, 2, global_mean)
  genoM0 <- sweep(geno$M, 2, global_mean)
  genoS0 <- sweep(geno$S, 2, global_mean)

  a <- t(genoN0) %*% V_invN %*% yN0 +
    t(genoM0) %*% V_invM %*% yM0 +
    t(genoS0) %*% V_invS %*% yS0
  c_val <- colSums(genoN0^2) + colSums(genoM0^2) + colSums(genoS0^2)

  set.seed(seed)
  n_snp <- ncol(genoN0)
  ratio_snp <- sample(seq_len(n_snp), min(n_ratio_snp, n_snp), replace = FALSE)
  genoN1 <- genoN0[, ratio_snp, drop = FALSE]
  genoM1 <- genoM0[, ratio_snp, drop = FALSE]
  genoS1 <- genoS0[, ratio_snp, drop = FALSE]

  b1 <- colSums(genoN1 * (V_invN %*% genoN1)) +
    colSums(genoM1 * (V_invM %*% genoM1)) +
    colSums(genoS1 * (V_invS %*% genoS1))
  c1 <- colSums(genoN1^2) + colSums(genoM1^2) + colSums(genoS1^2)
  gamma <- mean(b1 / c1)

  var <- 1 / (c_val * gamma)
  beta <- as.numeric(var * a)
  se <- sqrt(var)
  z <- beta / se

  if (is.null(snp)) {
    snp <- colnames(geno$M)
  }
  out <- data.frame(SNP = snp, beta = beta, var = as.numeric(var))
  out$pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  attr(out, "gamma") <- gamma
  out
}
