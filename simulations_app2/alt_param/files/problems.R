
dgp_null1 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- 0.4*Z[,1] + 0.5*Z[,2] + 0.6*Z[,3] - 0.7*Z[,4] + 1*Z[,7] + eps_X
  Y <- 0.6*Z[,1] - 0.2*Z[,2] + 0.3*Z[,4] + 0.9*Z[,5] - 0.5*Z[,6] + eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_null2 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- tanh(0.5*Z[,1] - 0.9*Z[,2] + Z[,3] + eps_X)
  Y <- exp(-0.8*Z[,4]*Z[,5] + 0.6*Z[,6]*Z[,7] + eps_Y)
  
  list(X = X, Y = Y, Z = Z)
}

dgp_null3 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- 0.5*Z[,1] - 0.9*Z[,2] + 0.4*Z[,3]^2 + (Z[,4]*Z[,5])*eps_X
  Y <- -0.8*Z[,4] + Z[,5]^2 + exp(Z[,6]) + sin(2*pi*Z[,7])*eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_null4 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- sin(2*pi*Z[,1]) + 0.1*eps_X
  Y <- sin(2*pi*Z[,1]) + eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_alt1 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- 0.7*Z[,1] + Z[,2] + eps_X
  Y <- 0.4*Z[,3] - 0.2*Z[,4] - 0.1*X + eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_alt2 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- sin(Z[,1]) + eps_X
  Y <- tanh(Z[,2]) + 0.2*X^2*Z[,3] + eps_Y
  
  list(X = X, Y = Y, Z = Z)
}

dgp_alt3 <- function(n, ...){
  
  eps_X <- rnorm(n, 0, 1)
  eps_Y <- rnorm(n, 0, 1)
  
  Z <- matrix(rnorm(n*7,0,1),ncol = 7)
  X <- 0.2*Z[,2]^3 + tanh(Z[,4]) + eps_X
  Y <- sin(pi*Z[,1]) - 0.4*Z[,2]^2 + cos(0.2*pi*X)*eps_Y
  
  list(X = X, Y = Y, Z = Z)
  
}
