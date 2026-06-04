compute_g_tilde_three_site <- function(gN, gM, gS,
                                       XXVX_invXVN, XXVX_invXVM, XXVX_invXVS,
                                       N_n, N_m, N_s) {
  gN <- as.matrix(gN)
  gM <- as.matrix(gM)
  gS <- as.matrix(gS)

  Zg <- XXVX_invXVN %*% gN + XXVX_invXVM %*% gM + XXVX_invXVS %*% gS
  Zg_N <- Zg[seq_len(N_n), , drop = FALSE]
  Zg_M <- Zg[(N_n + 1):(N_n + N_m), , drop = FALSE]
  Zg_S <- Zg[(N_n + N_m + 1):(N_n + N_m + N_s), , drop = FALSE]

  rbind(gN - Zg_N, gM - Zg_M, gS - Zg_S)
}

estimate_binary_grammar_gamma_ratio <- function(genoN, genoM, genoS,
                                                P_null, W_null,
                                                XXVX_invXVN, XXVX_invXVM, XXVX_invXVS,
                                                N_n, N_m, N_s,
                                                n_ratio_snp = 30,
                                                seed = 2025) {
  n_snp <- ncol(genoM)
  n_ratio_snp <- min(n_ratio_snp, n_snp)
  set.seed(seed)
  ratio_snp <- sample(seq_len(n_snp), n_ratio_snp, replace = FALSE)
  ratio_values <- numeric(n_ratio_snp)

  for (k in seq_along(ratio_snp)) {
    j <- ratio_snp[k]
    g_tilde <- compute_g_tilde_three_site(
      genoN[, j, drop = FALSE],
      genoM[, j, drop = FALSE],
      genoS[, j, drop = FALSE],
      XXVX_invXVN,
      XXVX_invXVM,
      XXVX_invXVS,
      N_n,
      N_m,
      N_s
    )
    exact_var <- as.numeric(t(g_tilde) %*% P_null %*% g_tilde)
    approx_var <- as.numeric(t(g_tilde) %*% (as.numeric(W_null) * g_tilde))
    ratio_values[k] <- exact_var / approx_var
  }

  ratio_values <- ratio_values[is.finite(ratio_values) & ratio_values > 0]
  if (length(ratio_values) == 0) {
    stop("No valid GRAMMAR-Gamma variance ratios were obtained.")
  }
  mean(ratio_values)
}

score_binary_three_site <- function(result, geno, snp = NULL, n_ratio_snp = 30, seed = 2025) {
  N_n <- length(result$site_y$N)
  N_m <- length(result$site_y$M)
  N_s <- length(result$site_y$S)

  res <- result$residuals
  resN <- as.matrix(res[seq_len(N_n)], ncol = 1)
  resM <- as.matrix(res[(N_n + 1):(N_n + N_m)], ncol = 1)
  resS <- as.matrix(res[(N_n + N_m + 1):(N_n + N_m + N_s)], ncol = 1)

  score <- suoscore(geno$N, geno$M, geno$S, resN, resM, resS)

  Z <- result$obj.noK$Z$XXVX_invXV
  V <- result$obj.noK$V
  XXVX_invXVN <- Z[, seq_len(N_n), drop = FALSE]
  XXVX_invXVM <- Z[, (N_n + 1):(N_n + N_m), drop = FALSE]
  XXVX_invXVS <- Z[, (N_n + N_m + 1):(N_n + N_m + N_s), drop = FALSE]

  var_ratio <- estimate_binary_grammar_gamma_ratio(
    genoN = geno$N,
    genoM = geno$M,
    genoS = geno$S,
    P_null = result$wen$suo$P,
    W_null = V,
    XXVX_invXVN = XXVX_invXVN,
    XXVX_invXVM = XXVX_invXVM,
    XXVX_invXVS = XXVX_invXVS,
    N_n = N_n,
    N_m = N_m,
    N_s = N_s,
    n_ratio_snp = n_ratio_snp,
    seed = seed
  )

  var <- parallel_fedgg(
    ncol(geno$M),
    geno$M,
    geno$N,
    geno$S,
    XXVX_invXVM,
    XXVX_invXVN,
    XXVX_invXVS,
    V,
    var_ratio,
    N_m,
    N_n,
    N_s
  )

  if (is.null(snp)) {
    snp <- colnames(geno$M)
  }
  out <- data.frame(SNP = snp, T = as.numeric(score), var = as.numeric(var))
  out$pval <- stats::pchisq((out$T)^2 / out$var, lower.tail = FALSE, df = 1)
  attr(out, "varRatio_null") <- var_ratio
  out
}
