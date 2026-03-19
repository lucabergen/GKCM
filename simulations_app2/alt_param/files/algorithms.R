
GKCM_RF_wrapper <- function(data, job, instance, ...){
  
  gaussian_kernel <- function(x, y = NULL, sigma = 1) {
    
    if(!is.vector(x)) x <- as.vector(x)
    
    if (is.null(y)) y <- x
    
    exp(-outer(x, y, "-")^2 / (2 * sigma^2))
    
  }
  
  center_kernel <- function(K) {
    
    n <- nrow(K)
    H <- diag(n) - matrix(1/n, n, n)
    
    H %*% K %*% H
    
  }
  
  compute_residuals <- function(K, W) {
    
    n <- nrow(K)
    
    # Compute B = (I - W)
    B <- diag(n) - W
    R <- B %*% K %*% t(B)
    
    # Symmetrize
    (R + t(R)) / 2
    
  }
  
  GKCM_RF <- function(X, Y, Z, seed = NULL, d = NULL, sigma = 1) {
    
    n <- nrow(Z)
    set.seed(seed)
    
    X <- scale(X)
    Y <- scale(Y)
    Z <- scale(Z)
    
    K <- gaussian_kernel(X, sigma = sigma)
    L <- gaussian_kernel(Y, sigma = sigma)
    
    # Set honesty = F to avoid sample splitting
    DRF_X <- drf::drf(Z, X, num.trees = ncol(Z)*100, bandwidth = sigma, 
                      honesty = F, response.scaling = F, mtry = 50,
                      min.node.size = 5)
    
    DRF_Y <- drf::drf(Z, Y, num.trees = ncol(Z)*100, bandwidth = sigma, 
                      honesty = F, response.scaling = F, mtry = 50, 
                      min.node.size = 5)
    
    W_X <- drf::get_sample_weights(DRF_X, newdata = Z) |> as.matrix()
    W_Y <- drf::get_sample_weights(DRF_Y, newdata = Z) |> as.matrix()
    
    K_XZ <- compute_residuals(K, W_X) |> center_kernel()
    L_YZ <- compute_residuals(L, W_Y) |> center_kernel()
    
    # Compute test statistic
    t_stat <- sum(K_XZ * L_YZ) / n
    
    G <- K_XZ * L_YZ
    H <- center_kernel(G)/(n - 1)
    
    ev <- eigen(H, symmetric = TRUE)$values
    ev <- ev[ev > 1e-16]
    
    if (is.null(d)){
      # Moment based approximation of the null  
      p_val <- tryCatch({
        1 - momentchi2::lpb4(ev, t_stat)
      }, error = function(e) {
        1 - momentchi2::hbe(ev, t_stat)
      })
    } else {
      # Sampling based approximation of the null
      m <- length(ev)
      Chisq <- matrix(rchisq(m * d, df = 1), nrow = m, ncol = d)
      # Elementwise multiply each column with EVs
      t_null <- colSums(Chisq * ev)  
      p_val <- mean(t_null >= t_stat)
    }
    
    list(t_stat = t_stat, p_val = p_val) 
    
  }
  
  GKCM_RF(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

GKCM_KRR_wrapper <- function(data, job, instance, ...){
  
  inv <- function(X) {chol2inv(chol(X))}
  
  gaussian_kernel <- function(X, sigma) {
    
    if (is.vector(X)) X <- matrix(X, ncol = 1)

    # Compute squared distances
    X2 <- rowSums(X^2)
    sq_dist <- outer(X2, X2, "+") - 2 * (X %*% t(X))
    # Censor
    sq_dist[sq_dist < 0] <- 0  
    
    exp(-sq_dist/(2 * sigma^2))
  }
  
  center_kernel <- function(K) {
    
    n <- nrow(K)
    H <- diag(n) - matrix(1/n, n, n)
    
    H %*% K %*% H
    
  }
  
  # L is the Gram matrix of the targets, K of the predictors
  loo_cv <- function(L, K, lambdas= c(1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1,1e+2,1e+3)) {
    
    n_full <- nrow(K)
    n <- n_full
    
    if (n_full > 1000) {
      idx <- sample.int(n_full, 1000)
      K <- K[idx, idx]
      L <- L[idx, idx]
      n <- 1000
    }
    
    I <- diag(n)
    
    # Symmetrize
    K <- 0.5 * (K + t(K))
    L <- 0.5 * (L + t(L))
    
    res <- numeric(length(lambdas))
    
    for (j in seq_along(lambdas)) {
      
      # Convert full-sample lambda to reduced-sample lambda
      lam <- lambdas[j] * (n / n_full)
      
      G <- K + lam * I 
      
      Ginv <- inv(G)
      
      A <- K %*% Ginv
      # Symmetrize
      A <- 0.5 * (A + t(A))
      
      E <- I - A
      
      num <- diag(E %*% L %*% E)
      den <- (1 - diag(A))^2
      
      res[j] <- sum(num / den)
    }
    
    list(lambda_opt = lambdas[which.min(res)], 
         loo_errors = setNames(res, lambdas))
  }
  
  GKCM_KRR <- function(X, Y, Z, lambda = NULL, sigma = 1, 
                       seed = NULL, d = NULL) {
    
    X <- scale(X)
    Y <- scale(Y)
    Z <- scale(Z)
    
    med_heur <- median(dist(Z))
    
    n <- nrow(Z)
    p <- ncol(Z)
    set.seed(seed)
    
    K <- gaussian_kernel(X, sigma = sigma)
    L <- gaussian_kernel(Y, sigma = sigma)
    M <- gaussian_kernel(Z, med_heur)
    
    if (!is.null(lambda)){
      lam_X <- lambda
      lam_Y <- lambda
    } else {
      lam_X <- loo_cv(K, M)$lambda_opt
      lam_Y <- loo_cv(L, M)$lambda_opt
    }
    
    R_X <- lam_X * inv(M + lam_X * diag(n))
    R_Y <- lam_Y * inv(M + lam_Y * diag(n))
    
    # (Centered) matrices of inner products of residuals
    K_XZ <- R_X %*% K %*% R_X |> center_kernel()
    L_YZ <- R_Y %*% L %*% R_Y |> center_kernel()
    
    # Compute test statistic
    t_stat <- sum(K_XZ * L_YZ) / n
    
    G <- K_XZ * L_YZ
    H <- center_kernel(G)/(n - 1)
    
    ev <- eigen(H, symmetric = TRUE)$values
    ev <- ev[ev > 1e-16]
    
    if (is.null(d)){
      # Moment based approximation of the null  
      p_val <- tryCatch({
        1 - momentchi2::lpb4(ev, t_stat)
      }, error = function(e) {
        1 - momentchi2::hbe(ev, t_stat)
      })
    } else {
      # Sampling based approximation of the null
      m <- length(ev)
      Chisq <- matrix(rchisq(m * d, df = 1), nrow = m, ncol = d)
      # Elementwise multiply each column with EVs
      t_null <- colSums(Chisq * ev)  
      p_val <- mean(t_null >= t_stat)
    }
    
    list(t_stat = t_stat, p_val = p_val)
  }
  
  GKCM_KRR(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

PCM_wrapper <- function(data, job, instance, ...){
  
  args <- list(num.trees = 700, mtry = 7, min.node.size = 5)
  
  comets::pcm(
    Y = instance$X,
    X = instance$Y,
    Z = instance$Z,
    args_YonXZ = args,
    args_YonZ = args,
    args_YhatonZ = args,
    args_VonXZ = args,
    args_RonZ = args
  )$p.value
  
}

wGCM_wrapper <- function(data, job, instance, ...){
  
  args <- list(num.trees = 700, mtry = 7, min.node.size = 5)
  
  comets::wgcm(
    Y = instance$X,
    X = instance$Y,
    Z = instance$Z,
    args_YonZ = args,
    args_XonZ = args,
    args_wfun = args
  )$p.value[[1]]
  
}

GCM_wrapper <- function(data, job, instance, ...){

  args <- list(num.trees = 700, mtry = 7, min.node.size = 5)
  
  comets::gcm(
    Y = instance$X,
    X = instance$Y,
    Z = instance$Z,
    args_YonZ = args,
    args_XonZ = args
  )$p.value
  
}

KCIT_wrapper <- function(data, job, instance, ...){
  
  gaussian_kernel <- function(X, sigma) {
    
    if (is.vector(X)) X <- matrix(X, ncol = 1)
    
    # Compute squared distances
    X2 <- rowSums(X^2)
    sq_dist <- outer(X2, X2, "+") - 2 * (X %*% t(X))
    # Censor
    sq_dist[sq_dist < 0] <- 0  
    
    exp(-sq_dist/(2 * sigma^2))
  }
  
  center_kernel <- function(K) {
    
    n <- nrow(K)
    H <- diag(n) - matrix(1/n, n, n)
    
    H %*% K %*% H
    
  }
  
  inv <- function(X) {chol2inv(chol(X))}
  
  # L is the Gram matrix of the targets, K of the predictors
  loo_cv <- function(L, K, lambdas= c(1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1,1e+2,1e+3)) {
    
    n_full <- nrow(K)
    n <- n_full
    
    if (n_full > 1000) {
      idx <- sample.int(n_full, 1000)
      K <- K[idx, idx]
      L <- L[idx, idx]
      n <- 1000
    }
    
    I <- diag(n)
    
    # Symmetrize
    K <- 0.5 * (K + t(K))
    L <- 0.5 * (L + t(L))
    
    res <- numeric(length(lambdas))
    
    for (j in seq_along(lambdas)) {
      
      # Convert full-sample lambda to reduced-sample lambda
      lam <- lambdas[j] * (n / n_full)
      
      G <- K + lam * I 
      
      Ginv <- inv(G)
      
      A <- K %*% Ginv
      # Symmetrize
      A <- 0.5 * (A + t(A))
      
      E <- I - A
      
      num <- diag(E %*% L %*% E)
      den <- (1 - diag(A))^2
      
      res[j] <- sum(num / den)
    }
    
    list(lambda_opt = lambdas[which.min(res)], 
         loo_errors = setNames(res, lambdas))
  }
  
  # U-statistic from He et al., appendix H.1
  u_stat <- function(K, L) {
    
    n <- nrow(K)
    
    # Set to zero to compute statistic only with i != j as in
    # https://github.com/he-zh/kci-hardness/blob/main/u_estimator.py
    diag(K) <- 0
    diag(L) <- 0
    
    o <- rep(1, n)
    
    # trace term
    t1 <- sum(K * L)
    
    t2 <- as.numeric(((t(o) %*% K %*% o)*(t(o) %*% L %*% o)) / ((n-1) * (n-2)))
    
    t3 <- as.numeric((2 * (t(o) %*% K %*% L %*% o)) / (n - 1))
    
    as.numeric((t1 + t2 - t3) / (n * (n - 3)))
  }
  
  KCIT <- function(X, Y, Z, lambda = NULL, sigma = 1, B = 500, seed = NULL) {
    
    X <- scale(X)
    Y <- scale(Y)
    Z <- scale(Z)

    n <- nrow(Z)
    
    med_heur <- median(dist(Z))
    
    set.seed(seed)
    
    K <- gaussian_kernel(X, sigma = sigma)
    L <- gaussian_kernel(Y, sigma = sigma)
    M <- gaussian_kernel(Z, sigma = med_heur)
    
    if (!is.null(lambda)){
      lam_X <- lambda
      lam_Y <- lambda
    } else {
      # LOO CV for lambda (on max. 1000 observation)
      lam_X <- loo_cv(K, M)$lambda_opt
      lam_Y <- loo_cv(L, M)$lambda_opt
    }
    
    R_X <- lam_X * inv(M + lam_X * diag(n))
    R_Y <- lam_Y * inv(M + lam_Y * diag(n))
    
    # Residual Gram matrices
    K_c <- R_X %*% K %*% R_X
    L_c <- R_Y %*% L %*% R_Y
    
    # Decomposition of joint embedding of (X,Z) via weighting; 
    # set higher bandwidth to minimize type-I error as described by He et al.
    L_c_w <- L_c * gaussian_kernel(Z, sigma = 2*med_heur)
    
    t_stat <- u_stat(K_c, L_c_w)
    
    # Wild bootstrap from He et al., appendix H.1
    t_null <- numeric(B)
    for (b in seq_len(B)) {
      q <- sample(c(-1, 1), size = n, replace = TRUE)
      K_c_null <- outer(q, q) * K_c
      t_null[b] <- u_stat(K_c_null, L_c_w)
    }
    
    p_val <- (1 + sum(t_null > t_stat)) / (B + 1)
    
    
    list(t_stat = t_stat, p_val = p_val)
  }
  
  KCIT(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

RCIT_wrapper <- function(data, job, instance, ...){
  
  RFF <- function(X, n_feat, sigma = NULL){ 
    
    if (is.vector(X)){X <- as.matrix(X, ncol = 1)} 
    
    n = nrow(X) 
    p = ncol(X) 
    
    if(length(sigma) < p) { 
      if (length(sigma) > 1) {
        stop("Length of sigma vector does not match number of columns") 
      } else {
        sigma <- rep(sigma, p)
      }
    }
    
    # Compute RFFs 
    W <- matrix(rnorm(n_feat * p), ncol = p) 
    W <- sweep(W, 2, sigma, "/") 
    
    b <- runif(n_feat, min = 0, max = 2*pi) 
    
    M <- sweep(X %*% t(W), 2, b, "+") 
    
    sqrt(2/n_feat) * cos(M) 
    
  } 
  
  inv <- function(X) {chol2inv(chol(X))}
  
  loo_cv <- function(Y, X, lambdas = c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2)) {
    
    n_full <- nrow(X)
    
    if (n_full > 1000) {
      idx <- sample.int(n_full, 1000)
      X <- X[idx, , drop = FALSE]
      Y <- Y[idx, , drop = FALSE]
    }
    
    n <- nrow(X)
    p <- ncol(X)
    Ip <- diag(p)
    
    XtX <- crossprod(X)
    XtY <- crossprod(X, Y)
    
    res <- numeric(length(lambdas))
    
    for (j in seq_along(lambdas)) {
      # Convert full-sample lambda to reduced-sample lambda
      lam <- lambdas[j] * (n / n_full)
      
      G <- XtX + lam * Ip
      
      Ginv <- tryCatch({
        inv(G)
      }, error = function(e) {
        solve(G)
      })
      
      # Compute residuals
      B <- Ginv %*% XtY
      Yhat <- X %*% B
      R <- Y - Yhat
      
      s <- rowSums((X %*% Ginv) * X)
      
      num <- rowSums(R^2)
      den <- (1 - s)^2
      
      res[j] <- sum(num / den)
    }
    
    list(
      lambda_opt  = lambdas[which.min(res)],
      loo_errors  = setNames(res, lambdas)
    )
  }
  
  RCIT <- function(X, Y, Z, lambda = NULL, sigma = 1,  
                   n_RFF_xy = 10, n_RFF_z = 200, seed = NULL){
    
    set.seed(seed)
    
    X <- scale(X)
    Y <- scale(Y)
    Z <- scale(Z)
    
    n <- nrow(Z)
    med_heur <- median(dist(Z))
    
    XZ_RFF <- RFF(cbind(X,Z), n_RFF_xy, 
                 sigma = c(sigma, rep(2*med_heur, ncol(Z))))
    Y_RFF <- RFF(Y, n_RFF_xy, sigma = sigma)
    Z_RFF <- RFF(Z, n_RFF_z, sigma = med_heur)
    
    XZ_RFF <- scale(XZ_RFF, scale = F)
    Y_RFF <- scale(Y_RFF, scale = F)
    Z_RFF <- scale(Z_RFF, scale = F)
    
    if (!is.null(lambda)){
      lam_X <- lambda
      lam_Y <- lambda
    } else {
      lam_X <- loo_cv(XZ_RFF, Z_RFF)$lambda_opt
      lam_Y <- loo_cv(Y_RFF, Z_RFF)$lambda_opt
    }
    
    res_x <- KRR_residuals(XZ_RFF, Z_RFF, lam_X)
    res_y <- KRR_residuals(Y_RFF,  Z_RFF, lam_Y)
    
    res_x <- scale(res_x, scale = FALSE)
    res_y <- scale(res_y, scale = FALSE)
    
    m <- ncol(res_x)
    q <- ncol(res_y)
    
    C <- cov(res_x, res_y)
    t_stat <- n * sum(C^2)
    
    res <- do.call(cbind, lapply(1:q, function(j) res_x * res_y[, j]))
    
    Cov <- crossprod(res)/n
    ev  <- eigen(Cov, symmetric = T, only.values = T)$values
    ev  <- ev[ev > 1e-12]
    
    p_val <- tryCatch({
      1 - momentchi2::lpb4(ev, t_stat)
    }, error = function(e) {
      1 - momentchi2::hbe(ev, t_stat)
    })
    
    list(t_stat = t_stat, p_val = p_val)
    
  }
  
  RCIT(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

RCoT_wrapper <- function(data, job, instance, ...){
  
  RFF <- function(X, n_feat, sigma = NULL){ 
    
    if (is.vector(X)){X <- as.matrix(X, ncol = 1)} 
    
    n = nrow(X) 
    p = ncol(X) 
    
    if(length(sigma) < p) { 
      if (length(sigma) > 1) {
        stop("Length of sigma vector does not match number of columns") 
      } else {
        sigma <- rep(sigma, p)
      }
    }
    
    # Compute RFFs 
    W <- matrix(rnorm(n_feat * p), ncol = p) 
    W <- sweep(W, 2, sigma, "/") 
    
    b <- runif(n_feat, min = 0, max = 2*pi) 
    
    M <- sweep(X %*% t(W), 2, b, "+") 
    
    sqrt(2/n_feat) * cos(M) 
    
  } 
  
  inv <- function(X) {chol2inv(chol(X))}
  
  loo_cv <- function(Y, X, lambdas = c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2)) {
    
    n_full <- nrow(X)
    
    if (n_full > 1000) {
      idx <- sample.int(n_full, 1000)
      X <- X[idx, , drop = FALSE]
      Y <- Y[idx, , drop = FALSE]
    }
    
    n <- nrow(X)
    p <- ncol(X)
    Ip <- diag(p)
    
    XtX <- crossprod(X)
    XtY <- crossprod(X, Y)
    
    res <- numeric(length(lambdas))
    
    for (j in seq_along(lambdas)) {
      # Convert full-sample lambda to reduced-sample lambda
      lam <- lambdas[j] * (n / n_full)
      
      G <- XtX + lam * Ip
      
      Ginv <- tryCatch({
        inv(G)
      }, error = function(e) {
        solve(G)
      })
      
      # Compute residuals
      B <- Ginv %*% XtY
      Yhat <- X %*% B
      R <- Y - Yhat
      
      s <- rowSums((X %*% Ginv) * X)
      
      num <- rowSums(R^2)
      den <- (1 - s)^2
      
      res[j] <- sum(num / den)
    }
    
    list(
      lambda_opt  = lambdas[which.min(res)],
      loo_errors  = setNames(res, lambdas)
    )
  }
  
  # Compute residuals
  KRR_residuals <- function(Y, X, lambda){
    
    n <- nrow(X)
    p <- ncol(X)
    
    XtX <- crossprod(X)
    diag(XtX) <- diag(XtX) + lambda
    R   <- chol(XtX)
    
    XtY <- crossprod(X, Y)
    B   <- backsolve(R, forwardsolve(t(R), XtY))
    
    Y - X %*% B
    
  }
  
  RCoT <- function(X, Y, Z, lambda = NULL, sigma = 1, 
                   n_RFF_xy = 10, n_RFF_z = 200, seed = NULL){
  
  set.seed(seed)
    
  X <- scale(X)
  Y <- scale(Y)
  Z <- scale(Z)
  
  n <- nrow(Z)
  med_heur <- median(dist(Z))
  
  X_RFF <- RFF(X, n_RFF_xy, sigma = sigma)
  Y_RFF <- RFF(Y, n_RFF_xy, sigma = sigma)
  Z_RFF <- RFF(Z, n_RFF_z, sigma = med_heur)
  
  X_RFF <- scale(X_RFF, scale = F)
  Y_RFF <- scale(Y_RFF, scale = F)
  Z_RFF <- scale(Z_RFF, scale = F)
  
  if (!is.null(lambda)){
    lam_X <- lambda
    lam_Y <- lambda
  } else {
    lam_X <- loo_cv(X_RFF, Z_RFF)$lambda_opt
    lam_Y <- loo_cv(Y_RFF, Z_RFF)$lambda_opt
  }
  
  res_x <- KRR_residuals(X_RFF, Z_RFF, lam_X)
  res_y <- KRR_residuals(Y_RFF, Z_RFF, lam_Y)
  
  res_x <- scale(res_x, scale = FALSE)
  res_y <- scale(res_y, scale = FALSE)
  
  m <- ncol(res_x)
  q <- ncol(res_y)
  
  C <- cov(res_x, res_y)
  t_stat <- n * sum(C^2)
  
  res <- do.call(cbind, lapply(1:q, function(j) res_x * res_y[, j]))
  
  Cov <- crossprod(res)/n
  ev  <- eigen(Cov, symmetric = T, only.values = T)$values
  ev  <- ev[ev > 1e-12]
  
  p_val <- tryCatch({
    1 - momentchi2::lpb4(ev, t_stat)
  }, error = function(e) {
    1 - momentchi2::hbe(ev, t_stat)
  })
  
  list(t_stat = t_stat, p_val = p_val)
  
  }
  
  RCoT(
    X = instance$X,
    Y = instance$Y,
    Z = instance$Z
  )$p_val
  
}

