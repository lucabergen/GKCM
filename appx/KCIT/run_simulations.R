
# Install necessary packages
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(c("future@1.69.0", "future.apply@1.20.1", "here@1.0.2", 
                   "data.table@1.18.2.1", "comets@0.2-2", "drf@1.1.0", 
                   "momentchi2@0.1.5"))

library(future)
library(future.apply)
library(here)
library(data.table)

source(here("appx/KCIT/CI_tests.R"))

# number of iterations
k <- 100  
# number of cores
n_cores <- max(1, parallel::detectCores() - 1)

alg_names <- c("GKCM_RF","GKCM_KRR","PCM","wGCM","GCM","KCIT","RCIT","RCoT")

future_globals <- c(
  "alg_names", "alg_fns", "dims", "eval_rep",
  "dgp_caseI_CI", "dgp_caseI_CD", "dgp_caseII_CI", "dgp_caseII_CD", "zscore_vec",
  "gaussian_kernel", "center_kernel", "compute_residuals",
  "u_stat", "RFF", "KRR_residuals",
  "GKCM_RF", "GKCM_KRR", "PCM", "wGCM", "GCM", "KCIT", "RCIT", "RCoT"
)

# Study parameters per Zhang et al.
set.seed(1)
n_vals <- c(200, 400)
dims   <- 1:5
alphas <- c(0.01, 0.05)

plan(multisession, workers = n_cores)

stopifnot(is.character(alg_names), length(alg_names) == 8)
alg_fns <- mget(alg_names, mode = "function", inherits = TRUE)

zscore_vec <- function(v) (v - mean(v)) / sd(v)

eval_rep <- function(dat_by_d) {
  pmat <- matrix(
    NA_real_,
    nrow = length(dims),
    ncol = length(alg_names),
    dimnames = list(d = as.character(dims), alg = alg_names)
  )
  
  for (d in dims) {
    x <- dat_by_d[[d]]$x
    y <- dat_by_d[[d]]$y
    Z <- dat_by_d[[d]]$z
    
    for (alg in alg_names) {
      pmat[as.character(d), alg] <- tryCatch(
        as.numeric(alg_fns[[alg]](x, y, Z)),
        error = function(e) NA_real_
      )
    }
  }
  pmat
}

pmats_to_df <- function(pmats, case, CI, n) {
  arr <- simplify2array(pmats)
  
  # enforce dimnames
  dimnames(arr) <- list(
    d    = as.character(dims),
    alg  = alg_names,
    iter = as.character(seq_along(pmats))
  )
  
  df <- as.data.table(as.table(arr))
  setnames(df, c("d", "alg_name", "iter", "p_value"))
  
  df[, `:=`(
    d    = as.integer(as.character(d)),
    iter = as.integer(as.character(iter)),
    case = case,
    CI   = CI,
    n    = n
  )]
  
  df[, .(case, CI, n, iter, d, alg_name, p_value)]
}
# Case I DGPs (one replication)
dgp_caseI_CI <- function(n) {
  x0 <- rnorm(n); y0 <- rnorm(n); z0 <- rnorm(n)
  
  zz1 <- 0.7 * ((z0^3) / 5 + z0 / 2)
  x <- zz1 + tanh(x0)
  x <- x + (x^3) / 3 + tanh(x / 3) / 2
  
  zz2 <- (z0^3 / 4 + z0) / 3
  y <- y0 + zz2
  y <- y + tanh(y / 3)
  
  x <- zscore_vec(x); y <- zscore_vec(y); z1 <- zscore_vec(z0)
  
  out <- vector("list", 5)
  for (d in 1:5) {
    Z <- if (d == 1) matrix(z1, ncol = 1) else cbind(z1, matrix(rnorm(n * (d - 1)), nrow = n))
    out[[d]] <- list(x = x, y = y, z = Z)
  }
  out
}

dgp_caseI_CD <- function(n) {
  x0 <- rnorm(n); y0 <- rnorm(n); z0 <- rnorm(n)
  
  zz1 <- 0.7 * ((z0^3) / 5 + z0 / 2)
  x <- zz1 + tanh(x0)
  x <- x + (x^3) / 3 + tanh(x / 3) / 2
  
  zz2 <- (z0^3 / 4 + z0) / 3
  y <- y0 + zz2
  y <- y + tanh(y / 3)
  
  x <- zscore_vec(x); y <- zscore_vec(y); z1 <- zscore_vec(z0)
  
  ff <- rnorm(n) * 0.5
  x <- x + ff
  y <- y + ff
  
  out <- vector("list", 5)
  for (d in 1:5) {
    Z <- if (d == 1) matrix(z1, ncol = 1) else cbind(z1, matrix(rnorm(n * (d - 1)), nrow = n))
    out[[d]] <- list(x = x, y = y, z = Z)
  }
  out
}

# Case II DGPs (one replication)
dgp_caseII_CI <- function(n) {
  x0 <- rnorm(n); y0 <- rnorm(n); z0 <- rnorm(n)
  
  zz1_raw <- 0.7 * ((z0^3) / 5 + z0 / 2)
  x1 <- zz1_raw + tanh(x0)
  x1 <- x1 + (x1^3) / 3 + tanh(x1 / 3) / 2
  
  zz2_raw <- (z0^3 / 4 + z0) / 3
  y1 <- y0 + zz2_raw
  y1 <- y1 + tanh(y1 / 3)
  
  x1 <- zscore_vec(x1); y1 <- zscore_vec(y1); z1 <- zscore_vec(z0)
  
  out <- vector("list", 5)
  out[[1]] <- list(x = x1, y = y1, z = matrix(z1, ncol = 1))
  
  phi <- function(t) t / 2 + 0.7 * tanh(t)
  
  Zmat <- matrix(z1, ncol = 1)
  zz1_prev <- NULL
  zz2_prev <- NULL
  
  for (d in 2:5) {
    zd_raw <- rnorm(n)
    x0 <- rnorm(n); y0 <- rnorm(n)
    
    if (d == 2) {
      base1 <- zz1_raw / 2 + zd_raw
      base2 <- zz2_raw / 2 + zd_raw
    } else {
      base1 <- zz1_prev * (2/3) + zd_raw * (5/6)
      base2 <- zz2_prev * (2/3) + zd_raw * (5/6)
    }
    
    zz1_d <- phi(base1)
    zz2_d <- phi(base2)
    
    x <- zz1_d + tanh(x0)
    x <- x + (x^3) / 3 + tanh(x / 3) / 2
    
    y <- y0 + zz2_d
    y <- y + tanh(y / 3)
    
    x <- zscore_vec(x); y <- zscore_vec(y)
    zd <- zscore_vec(zd_raw)
    Zmat <- cbind(Zmat, zd)
    
    out[[d]] <- list(x = x, y = y, z = Zmat)
    
    zz1_prev <- zz1_d
    zz2_prev <- zz2_d
  }
  
  out
}

dgp_caseII_CD <- function(n) {
  x0 <- rnorm(n); y0 <- rnorm(n); z0 <- rnorm(n)
  
  zz1_raw <- 0.7 * ((z0^3) / 5 + z0 / 2)
  x1 <- zz1_raw + tanh(x0)
  x1 <- x1 + (x1^3) / 3 + tanh(x1 / 3) / 2
  
  zz2_raw <- (z0^3 / 4 + z0) / 3
  y1 <- y0 + zz2_raw
  y1 <- y1 + tanh(y1 / 3)
  
  x1 <- zscore_vec(x1); y1 <- zscore_vec(y1); z1 <- zscore_vec(z0)
  
  ff <- rnorm(n) * 0.5
  x1 <- x1 + ff
  y1 <- y1 + ff
  
  out <- vector("list", 5)
  out[[1]] <- list(x = x1, y = y1, z = matrix(z1, ncol = 1))
  
  phi <- function(t) t / 2 + 0.7 * tanh(t)
  
  Zmat <- matrix(z1, ncol = 1)
  zz1_prev <- NULL
  zz2_prev <- NULL
  
  for (d in 2:5) {
    zd_raw <- rnorm(n)
    x0 <- rnorm(n); y0 <- rnorm(n)
    
    if (d == 2) {
      base1 <- zz1_raw / 2 + zd_raw
      base2 <- zz2_raw / 2 + zd_raw
    } else {
      base1 <- zz1_prev * (2/3) + zd_raw * (5/6)
      base2 <- zz2_prev * (2/3) + zd_raw * (5/6)
    }
    
    zz1_d <- phi(base1)
    zz2_d <- phi(base2)
    
    x <- zz1_d + tanh(x0)
    x <- x + (x^3) / 3 + tanh(x / 3) / 2
    
    y <- y0 + zz2_d
    y <- y + tanh(y / 3)
    
    x <- zscore_vec(x); y <- zscore_vec(y)
    zd <- zscore_vec(zd_raw)
    Zmat <- cbind(Zmat, zd)
    
    x <- x + ff
    y <- y + ff
    
    out[[d]] <- list(x = x, y = y, z = Zmat)
    
    zz1_prev <- zz1_d
    zz2_prev <- zz2_d
  }
  
  out
}

# Run scenarios
simulate_caseI_CI <- function(n, k) {
  pmats <- future_lapply(
    seq_len(k), 
    function(iter) eval_rep(dgp_caseI_CI(n)), 
    future.seed = TRUE, 
    future.globals = future_globals
  )
  pmats_to_df(pmats, case = "I", CI = TRUE, n = n)
}

simulate_caseI_CD <- function(n, k) {
  pmats <- future_lapply(
    seq_len(k), 
    function(iter) eval_rep(dgp_caseI_CD(n)), 
    future.seed = TRUE, 
    future.globals = future_globals
  )
  pmats_to_df(pmats, case = "I", CI = FALSE, n = n
  )
}

simulate_caseII_CI <- function(n, k) {
  pmats <- future_lapply(
    seq_len(k), function(iter) eval_rep(dgp_caseII_CI(n)), 
    future.seed = TRUE, 
    future.globals = future_globals
  )
  pmats_to_df(pmats, case = "II", CI = TRUE, n = n)
}

simulate_caseII_CD <- function(n, k) {
  pmats <- future_lapply(
    seq_len(k), 
    function(iter) eval_rep(dgp_caseII_CD(n)), 
    future.seed = TRUE, 
    future.globals = future_globals
  )
  pmats_to_df(pmats, case = "II", CI = FALSE, n = n)
}

# Run everything for n = 200, 400
res_I_CI  <- do.call(rbind, lapply(n_vals, simulate_caseI_CI,  k = k))
res_I_CD  <- do.call(rbind, lapply(n_vals, simulate_caseI_CD,  k = k))
res_II_CI <- do.call(rbind, lapply(n_vals, simulate_caseII_CI, k = k))
res_II_CD <- do.call(rbind, lapply(n_vals, simulate_caseII_CD, k = k))

res_all <- rbind(res_I_CI, res_I_CD, res_II_CI, res_II_CD)

save(res_all, file = here("appx/KCIT/results_df.Rdata"))
