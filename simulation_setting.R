library(MASS)
library(sparsepca)

################################################################################
########## Generate data for scenario 1~5 ######################################
###############################################################################

simulatedata <- function(seed, scenario, p = 50, noise_sd=0.1) {
  print(paste0("scenario ", scenario))
  set.seed(seed)
  n <- 100
  if(scenario %in% c(1, 2)) {
    V1 <- rnorm(n)
    V2 <- runif(n)
    X1 <- outer(V1, rep(1, 4)) + matrix(rnorm(n * 4, 0, 0.01), nrow = n)
    X2 <- outer(V2, rep(1, 2)) + matrix(rnorm(n * 2, 0, 0.05), nrow = n)
    X3 <- matrix(rnorm(n * (p - 6), 0, 0.0), nrow = n)
    XX <- cbind(scale(cbind(X1, X2)), X3)
    if(scenario %in% c(1)) {Xnoise <- 0.1}
    if(scenario %in% c(2)) {Xnoise <- 0.5}
    X <- XX + matrix(rnorm(n * p, 0, Xnoise), nrow = n)
    apply(X, 2, var)[1:6]
    spca(X, k = 2, alpha = 0.01, verbose = FALSE)$loadings[1:6, ]
  }
  if(scenario %in% c(3)) {
    V1 <- rnorm(n)
    V2 <- runif(n)
    V3 <- rnorm(n)
    X1 <- outer(V1, rep(1, 4)) + matrix(rnorm(n * 4, 0, 0.1), nrow = n)
    X2 <- outer(V2, rep(1, 2)) + matrix(rnorm(n * 2, 0, 0.1), nrow = n)
    X3 <- outer(V3, rep(1, 9)) + matrix(rnorm(n * 9, 0, 0.1), nrow = n)
    X4 <- matrix(rnorm(n * (p - 15), 0, 0.0), nrow = n)
    XX <- cbind(X1, X2, X3, X4)
    Xnoise <- 0.1  # 0.1, 0.3, 0.5
    X <- XX + matrix(rnorm(n * p, 0, Xnoise), nrow = n)
    apply(X, 2, var)[1:14]
    spca(X, k = 3, alpha = 0.01, verbose = FALSE)$loadings[1:15, ]
  }
  if(scenario %in% c(4)) {
    V1 <- rnorm(n)
    V2 <- runif(n)
    X1 <- outer(V1, rep(1, 4)) + matrix(rnorm(n * 4, 0, 0.1), nrow = n)
    X2 <- outer(V2, rep(1, 2)) + matrix(rnorm(n * 2, 0, 0.1), nrow = n)
    XX <- cbind(X1, X2)
    VV <- matrix(rnorm(n * 10), ncol = 10)
    for(j in 1:ncol(VV)) {
      XX <- cbind(XX, outer(VV[, j], rep(1, 4)) + matrix(rnorm(n * 4, 0, 0.1), nrow = n))
    }
    X3 <- matrix(rnorm(n * (p - ncol(XX)), 0, 0.0), nrow = n)
    XX <- cbind(XX, X3)
    Xnoise <- 0.1
    X <- XX + matrix(rnorm(n * p, 0, Xnoise), nrow = n)
    spca(X, k = 12, alpha = 0.01, verbose = FALSE)$loadings[, ]
  }
  if(scenario %in% c(5)) {
    rho <- 0.8; Sigma <- rho^abs(outer(1:p, 1:p, "-")) # AR(1)
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  }
  V <- matrix(c(0.5, 0.5, 0.5, 0.5, rep(0, p - 4), 
                0, 0, 0, 0, sqrt(0.5), sqrt(0.5), rep(0, p - 6)), nrow = p)
  Ucenters <- matrix(c(2 / sqrt(3), -4 / sqrt(3), 2 / sqrt(3), -2, 0, 2), ncol = 2) # different for every individual
  U <- Ucenters[sample(1:nrow(Ucenters), size = n, replace = TRUE), ]
  B <- U %*% t(V)
  y <- rowSums(X * B) + rnorm(n, 0, noise_sd)
  return(list(y = y, X = X, U = U, V = V, B = B))
}

