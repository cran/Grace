# This function calculates Grace coefficient estimates
# Author: Sen Zhao
# Email: sendavid7@gmail.com

grace <- function(Y, X, L, lambda.L, lambda.1 = 0, lambda.2 = 0, normalize.L = FALSE, K = 10){
  lambda.L <- unique(sort(lambda.L, decreasing = TRUE))
  lambda.1 <- unique(sort(lambda.1, decreasing = TRUE))
  lambda.2 <- unique(sort(lambda.2, decreasing = TRUE))

  ori.Y <- Y
  ori.X <- X
  if(!is.null(ncol(Y))){
    stop("Error: Y is not a vector.")
  }
  if(length(Y) != nrow(X)){
    stop("Error: Dimensions of X and Y do not match.")
  }
  if(!isSymmetric(L)){
    stop("Error: L is not a symmetric matrix.")
  }
  if(ncol(X) != ncol(L)){
    stop("Error: Dimensions of X and L do not match.")
  }
  if(min(lambda.L) < 0 | min(lambda.2) < 0 | min(lambda.1) < 0){
    stop("Error: Grace tuning parameters must be non-negative.")
  }
  if(min(lambda.L) == 0 & min(lambda.2) == 0){
    stop("Error: At least one of the grace tuning parameters must be positive.")
  }
  
  Y <- Y - mean(Y)  # Center Y
  n <- nrow(X)
  p <- ncol(X)
  scale.fac <- attr(scale(X), "scaled:scale")
  X <- scale(X)     # Standardize X
  
  if(normalize.L){
    diag(L)[diag(L) == 0] <- 1
    L <- diag(1 / sqrt(diag(L))) %*% L %*% diag(1 / sqrt(diag(L)))  # Normalize L
  }
  
  emin <- min(eigen(min(lambda.L) * L + min(lambda.2) * diag(p))$values)  # Minimum eigenvalue of the penalty weight matrix
  if(emin < 1e-5){
    stop("Error: The penalty matrix (lambda.L * L + lambda.2 * I) is not always positive definite for all tuning parameters. Consider increase the value of lambda.2.")
  }
  

  # If more than one tuning parameter is provided, perform K-fold cross-validation  
  if((length(lambda.L) > 1) | (length(lambda.1) > 1) | (length(lambda.2) > 1)){
    tun <- cvGrace(X, Y, L, lambda.L, lambda.1, lambda.2, K)
    lambda.L <- tun[1]
    lambda.1 <- tun[2]
    lambda.2 <- tun[3]  
  }
  
  # See Li & Li (2008) for reference
  Lnew <- lambda.L * L + lambda.2 * diag(p)
  eL <- eigen(Lnew)
  S <- eL$vectors %*% sqrt(diag(eL$values))
  l2star <- 1
  l1star <- lambda.1
  Xstar <- rbind(X, sqrt(l2star) * t(S)) / sqrt(1 + l2star)
  Ystar <- c(Y, rep(0, p))
  gammastar <- l1star / sqrt(1 + l2star) / 2 / (n + p)
  betahatstar <- glmnet(Xstar, Ystar, lambda = gammastar, intercept = FALSE, standardize = FALSE, thresh = 1e-11)$beta[, 1]
  betahat <- betahatstar / sqrt(1 + l2star)
  
  truebetahat <- betahat / scale.fac  # Scale back coefficient estimate
  truealphahat <- mean(ori.Y - ori.X %*% truebetahat)
  return(list(intercept = truealphahat, beta = truebetahat))
}
