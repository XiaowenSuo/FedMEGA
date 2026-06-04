compute_logistic_local <- function(X_local, y_local, alpha) {
  alpha <- as.numeric(alpha)
  eta <- as.vector(X_local %*% alpha)
  p <- pmin(pmax(stats::plogis(eta), 1e-8), 1 - 1e-8)
  grad <- as.numeric(crossprod(X_local, p - y_local))
  W <- p * (1 - p)
  H <- crossprod(X_local, X_local * W)
  list(grad = grad, Hessian = H)
}

compute_linear_local <- function(X_local, y_local, beta) {
  beta <- as.numeric(beta)
  resid <- as.vector(X_local %*% beta - y_local)
  grad <- as.numeric(2 * crossprod(X_local, resid))
  H <- 2 * crossprod(X_local)
  list(grad = grad, Hessian = H)
}

newton_update <- function(theta, grad, H, ridge = 1e-6) {
  theta <- as.numeric(theta)
  grad <- as.numeric(grad)
  H_reg <- H + diag(ridge, length(theta))
  step <- tryCatch(
    solve(H_reg, grad),
    error = function(e) qr.solve(H_reg, grad)
  )
  as.numeric(theta - step)
}

#' Federated null logistic regression for three sites
#'
#' @param XN,XM,XS Local covariate matrices for sites N, M, and S.
#' @param yN,yM,yS Local binary phenotypes coded as 0/1.
#' @param alpha_init Initial fixed-effect vector.
#' @param threshold Convergence threshold.
#' @param max_iter Maximum iterations.
#'
#' @return Fixed-effect estimates.
#' @export
fed_null_logistic_three_site <- function(XN, yN, XM, yM, XS, yS,
                                         alpha_init = rep(0, ncol(XN)),
                                         threshold = 1e-6,
                                         max_iter = 500) {
  alpha <- as.numeric(alpha_init)
  for (iter in seq_len(max_iter)) {
    N_local <- compute_logistic_local(XN, yN, alpha)
    M_local <- compute_logistic_local(XM, yM, alpha)
    S_local <- compute_logistic_local(XS, yS, alpha)

    global_grad <- N_local$grad + M_local$grad + S_local$grad
    global_Hessian <- N_local$Hessian + M_local$Hessian + S_local$Hessian
    alpha_new <- newton_update(alpha, global_grad, global_Hessian)

    if (max(abs(alpha_new - alpha)) < threshold) {
      return(alpha_new)
    }
    alpha <- alpha_new
  }
  warning("Maximum number of iterations reached without convergence.")
  alpha
}

#' Federated null linear regression for three sites
#'
#' @param XN,XM,XS Local covariate matrices for sites N, M, and S.
#' @param yN,yM,yS Local continuous phenotypes.
#' @param beta_init Initial fixed-effect vector.
#' @param threshold Convergence threshold.
#' @param max_iter Maximum iterations.
#'
#' @return Fixed-effect estimates.
#' @export
fed_null_linear_three_site <- function(XN, yN, XM, yM, XS, yS,
                                       beta_init = rep(0, ncol(XN)),
                                       threshold = 1e-6,
                                       max_iter = 500) {
  beta <- as.numeric(beta_init)
  for (iter in seq_len(max_iter)) {
    N_local <- compute_linear_local(XN, yN, beta)
    M_local <- compute_linear_local(XM, yM, beta)
    S_local <- compute_linear_local(XS, yS, beta)

    global_grad <- N_local$grad + M_local$grad + S_local$grad
    global_Hessian <- N_local$Hessian + M_local$Hessian + S_local$Hessian
    beta_new <- newton_update(beta, global_grad, global_Hessian)

    if (max(abs(beta_new - beta)) < threshold) {
      return(beta_new)
    }
    beta <- beta_new
  }
  warning("Maximum number of iterations reached without convergence.")
  beta
}
