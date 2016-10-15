# This function calculates Grace coefficients and p-values
# Author: Sen Zhao
# Email: sen-zhao@sen-zhao.com

grace.test <- function(Y, X, L = NULL, lambda.L = NULL, lambda.2 = 0, normalize.L = FALSE, eta = 0.05, C = 4 * sqrt(3), K = 10){
  if(is.null(L)){
    L <- matrix(0, nrow = p, ncol = p)
    lambda.L <- 0
  }
  lambda.L <- unique(sort(lambda.L, decreasing = TRUE))
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
  if(min(lambda.L) < 0 | min(lambda.2) < 0){
    stop("Error: Tuning parameters must be non-negative.")
  }
  if(min(lambda.L) == 0 & min(lambda.2) == 0 & length(lambda.L) == 1 & length(lambda.2) == 1){
    stop("Error: At least one of the tuning parameters must be positive.")
  }
  
  Y <- Y - mean(Y)  # Center Y
  n <- nrow(X)
  p <- ncol(X)
  scale.fac <- attr(scale(X), "scaled:scale")
  X <- scale(X)     # Standardize X
  
  if(normalize.L & !is.null(L)){
    diag(L)[diag(L) == 0] <- 1
    L <- diag(1 / sqrt(diag(L))) %*% L %*% diag(1 / sqrt(diag(L)))  # Normalize L
  }
  
  emin <- min(eigen(min(lambda.L) * L + min(lambda.2) * diag(p))$values)  # Minimum eigenvalue of the penalty weight matrix
  if(emin < 1e-5){
    stop("Error: The penalty matrix (lambda.L * L + lambda.2 * I) is not always positive definite for all tuning parameters. Consider increase the value of lambda.2.")
  }
  
  # If more than one tuning parameter is provided, perform K-fold cross-validation  
  if((length(lambda.L) > 1) | (length(lambda.2) > 1)){
    tun <- cvGrace(X, Y, L, lambda.L, 0, lambda.2, K = K)
    lambda.L <- tun[1]
    lambda.2 <- tun[3]
    print(paste("Tuning parameters selected by ", K, "-fold cross-validation:", sep = ""))
    print(paste("lambda.L = ", lambda.L, sep = ""))
    print(paste("lambda.2 = ", lambda.2, sep = ""))
  }
  
  betahat <- c(solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% t(X) %*% Y)   # Grace coefficient estimate
  truebetahat <- betahat / scale.fac  # Scale back coefficient estimate
  truealphahat <- mean(ori.Y - ori.X %*% truebetahat)
  
  sig.L <- scalreg(X, Y)$hsigma   # Error standard deviation from the scaled lasso
  lam <- sig.L * C * sqrt(log(p) / n) / 2 # Lasso initial tuning parameter
  beta.init <- glmnet(X, Y, lambda = lam, intercept = FALSE)$beta[, 1]   # Initial estimator
  se <- sig.L * sqrt(diag(solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% t(X) %*% X %*% solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)))) # Standard error
  bias <- solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% (lambda.L * L + lambda.2 * diag(p)) %*% beta.init  # Bias correction
  targ <- abs(solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% (lambda.L * L + lambda.2 * diag(p)))
  diag(targ) <- 0
  correct <- apply(targ, 1, max) * (log(p) / n)^(0.5 - eta)  # Gamma
  
  teststat <- betahat + bias   # Bias corrected test statistic
  Tstat <- (abs(teststat) - correct) * (abs(teststat) > correct)
  pval <- 2 * (1 - pnorm(abs(Tstat / se)))
  
  return(list(intercept = truealphahat, beta = truebetahat, pvalue = pval))
}
