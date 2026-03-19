
# Utils
gaussian_kernel <- function(X, sigma) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  X2 <- rowSums(X^2)
  sq_dist <- outer(X2, X2, "+") - 2 * (X %*% t(X))
  sq_dist[sq_dist < 0] <- 0
  exp(-sq_dist / (2 * sigma^2))
}

center_kernel <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1/n, n, n)
  H %*% K %*% H
}

compute_residuals <- function(K, W) {
  n <- nrow(K)
  B <- diag(n) - W
  R <- B %*% K %*% t(B)
  (R + t(R)) / 2
}

# U-statistic from He et al., appendix H.1
u_stat <- function(K, L) {
  
  n <- nrow(K)
  diag(K) <- 0
  diag(L) <- 0
  o <- rep(1, n)
  
  t1 <- sum(K * L)
  t2 <- as.numeric(((t(o) %*% K %*% o) * (t(o) %*% L %*% o)) / ((n - 1) * (n - 2)))
  t3 <- as.numeric((2 * (t(o) %*% K %*% L %*% o)) / (n - 1))
  
  as.numeric((t1 + t2 - t3) / (n * (n - 3)))
}

# Random Fourier Features
RFF <- function(X, n_feat, sigma = NULL) {
  
  if (is.vector(X)) X <- as.matrix(X, ncol = 1)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(sigma) < p) {
    if (length(sigma) > 1) stop("Length of sigma vector does not match number of columns")
    sigma <- rep(sigma, p)
  }
  
  W <- matrix(rnorm(n_feat * p), ncol = p)
  W <- sweep(W, 2, sigma, "/")
  
  b <- runif(n_feat, min = 0, max = 2*pi)
  M <- sweep(X %*% t(W), 2, b, "+")
  sqrt(2 / n_feat) * cos(M)
}

# Ridge regression residuals
KRR_residuals <- function(Y, X, lambda) {
  
  XtX <- crossprod(X)
  diag(XtX) <- diag(XtX) + lambda
  R <- chol(XtX)
  
  XtY <- crossprod(X, Y)
  B <- backsolve(R, forwardsolve(t(R), XtY))
  
  Y - X %*% B
}


# Methods 

GKCM_RF <- function(X, Y, Z, seed = NULL, d = NULL, sigma = 1) {

  n <- nrow(Z)
  if (!is.null(seed)) set.seed(seed)
  
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  K <- gaussian_kernel(X, sigma = sigma)
  L <- gaussian_kernel(Y, sigma = sigma)
  
  DRF_X <- drf::drf(
    Z, X,
    num.trees = ncol(Z) * 100,
    bandwidth = sigma,
    honesty = FALSE,
    response.scaling = FALSE,
    min.node.size = 5,
    mtry = 50
  )
  
  DRF_Y <- drf::drf(
    Z, Y,
    num.trees = ncol(Z) * 100,
    bandwidth = sigma,
    honesty = FALSE,
    response.scaling = FALSE,
    min.node.size = 5,
    mtry = 50
  )
  
  W_X <- as.matrix(drf::get_sample_weights(DRF_X, newdata = Z))
  W_Y <- as.matrix(drf::get_sample_weights(DRF_Y, newdata = Z))
  
  K_XZ <- center_kernel(compute_residuals(K, W_X))
  L_YZ <- center_kernel(compute_residuals(L, W_Y))
  
  t_stat <- sum(K_XZ * L_YZ) / n
  
  G <- K_XZ * L_YZ
  H <- center_kernel(G) / (n - 1)
  
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 1e-16]
  
  if (is.null(d)) {
    p_val <- tryCatch(
      1 - momentchi2::lpb4(ev, t_stat),
      error = function(e) 1 - momentchi2::hbe(ev, t_stat)
    )
  } else {
    m <- length(ev)
    Chisq <- matrix(rchisq(m * d, df = 1), nrow = m, ncol = d)
    t_null <- colSums(Chisq * ev)
    p_val <- mean(t_null >= t_stat)
  }
  
  as.numeric(p_val)
}

GKCM_KRR <- function(X, Y, Z, lambda = 0.001, sigma = 1, seed = NULL, d = NULL) {

  if (!is.null(seed)) set.seed(seed)
  
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  n <- nrow(Z)
  med_heur <- median(dist(Z))
  
  K <- gaussian_kernel(X, sigma = sigma)
  L <- gaussian_kernel(Y, sigma = sigma)
  
  M_X <- gaussian_kernel(Z, sigma = med_heur)
  M_Y <- gaussian_kernel(Z, sigma = med_heur)
  
  R_X <- lambda * solve(M_X + lambda * diag(n))
  R_Y <- lambda * solve(M_Y + lambda * diag(n))
  
  K_XZ <- center_kernel(R_X %*% K %*% R_X)
  L_YZ <- center_kernel(R_Y %*% L %*% R_Y)
  
  t_stat <- sum(K_XZ * L_YZ) / n
  
  G <- K_XZ * L_YZ
  H <- center_kernel(G) / (n - 1)
  
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 1e-16]
  
  if (is.null(d)) {
    p_val <- tryCatch(
      1 - momentchi2::lpb4(ev, t_stat),
      error = function(e) 1 - momentchi2::hbe(ev, t_stat)
    )
  } else {
    m <- length(ev)
    Chisq <- matrix(rchisq(m * d, df = 1), nrow = m, ncol = d)
    t_null <- colSums(Chisq * ev)
    p_val <- mean(t_null >= t_stat)
  }
  
  as.numeric(p_val)
}

PCM <- function(X, Y, Z) {
  comets::pcm(
    Y = X,
    X = Y,
    Z = Z,
    reg_YonXZ   = "tuned_rf",
    reg_YonZ    = "tuned_rf",
    reg_YhatonZ = "tuned_rf",
    reg_VonXZ   = "tuned_rf",
    reg_RonZ    = "tuned_rf"
  )$p.value
}

wGCM <- function(X, Y, Z) {
  comets::wgcm(
    Y = X,
    X = Y,
    Z = Z,
    reg_YonZ = "tuned_rf",
    reg_XonZ = "tuned_rf",
    reg_wfun = "tuned_rf"
  )$p.value[[1]]
}

GCM <- function(X, Y, Z) {
  comets::gcm(
    Y = X,
    X = Y,
    Z = Z,
    reg_YonZ = "tuned_rf",
    reg_XonZ = "tuned_rf"
  )$p.value
}

KCIT <- function(X, Y, Z, lambda = 0.001, sigma = 1, B = 500, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  n <- nrow(Z)
  med_heur <- median(dist(Z))
  
  K <- gaussian_kernel(X, sigma = sigma)
  L <- gaussian_kernel(Y, sigma = sigma)
  M <- gaussian_kernel(Z, sigma = med_heur)
  
  A <- M + lambda * diag(n)
  R <- lambda * chol2inv(chol(A))
  
  # Residual Gram matrices
  K_c <- R %*% K %*% R
  L_c <- R %*% L %*% R
  
  # Weighting kernel for Z
  L_c_w <- L_c * gaussian_kernel(Z, sigma = 2 * med_heur)
  
  t_stat <- u_stat(K_c, L_c_w)
  
  # Wild bootstrap
  t_null <- numeric(B)
  for (b in seq_len(B)) {
    q <- sample(c(-1, 1), size = n, replace = TRUE)
    K_c_null <- outer(q, q) * K_c
    t_null[b] <- u_stat(K_c_null, L_c_w)
  }
  
  p_val <- (1 + sum(t_null > t_stat)) / (B + 1)
  as.numeric(p_val)
}

RCIT <- function(X, Y, Z, lambda = 0.001, sigma = 1, n_RFF_xy = 10, n_RFF_z = 200, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  n <- nrow(Z)
  med_heur <- median(dist(Z))
  
  XZ_RFF <- scale(RFF(cbind(X, Z), n_RFF_xy, sigma = c(sigma, rep(2 * med_heur, ncol(Z)))))
  Y_RFF  <- scale(RFF(Y, n_RFF_xy, sigma = sigma))
  Z_RFF  <- scale(RFF(Z, n_RFF_z, sigma = med_heur))
  
  res_x <- KRR_residuals(XZ_RFF, Z_RFF, lambda)
  res_y <- KRR_residuals(Y_RFF,  Z_RFF, lambda)
  
  res_x <- scale(res_x, scale = FALSE)
  res_y <- scale(res_y, scale = FALSE)
  
  C <- cov(res_x, res_y)
  t_stat <- n * sum(C^2)
  
  q <- ncol(res_y)
  res <- do.call(cbind, lapply(seq_len(q), function(j) res_x * res_y[, j]))
  
  Cov <- crossprod(res) / n
  ev  <- eigen(Cov, symmetric = TRUE, only.values = TRUE)$values
  ev  <- ev[ev > 1e-12]
  
  p_val <- tryCatch(
    1 - momentchi2::lpb4(ev, t_stat),
    error = function(e) 1 - momentchi2::hbe(ev, t_stat)
  )
  
  as.numeric(p_val)
}

RCoT <- function(X, Y, Z, lambda = 0.001, sigma = 1, n_RFF_xy = 10, n_RFF_z = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  n <- nrow(Z)
  med_heur <- median(dist(Z))
  
  X_RFF <- scale(RFF(X, n_RFF_xy, sigma = sigma))
  Y_RFF <- scale(RFF(Y, n_RFF_xy, sigma = sigma))
  Z_RFF <- scale(RFF(Z, n_RFF_z, sigma = med_heur))
  
  res_x <- KRR_residuals(X_RFF, Z_RFF, lambda)
  res_y <- KRR_residuals(Y_RFF, Z_RFF, lambda)
  
  res_x <- scale(res_x, scale = FALSE)
  res_y <- scale(res_y, scale = FALSE)
  
  C <- cov(res_x, res_y)
  t_stat <- n * sum(C^2)
  
  q <- ncol(res_y)
  res <- do.call(cbind, lapply(seq_len(q), function(j) res_x * res_y[, j]))
  
  Cov <- crossprod(res) / n
  ev  <- eigen(Cov, symmetric = TRUE, only.values = TRUE)$values
  ev  <- ev[ev > 1e-12]
  
  p_val <- tryCatch(
    1 - momentchi2::lpb4(ev, t_stat),
    error = function(e) 1 - momentchi2::hbe(ev, t_stat)
  )
  
  as.numeric(p_val)
}
