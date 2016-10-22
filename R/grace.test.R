# This function calculates Grace coefficients and p-values
# Author: Sen Zhao
# Email: sen-zhao@sen-zhao.com
# ----------------------------------------------------------------------------
# Arguments:
# Y: n by 1 vector of the response variable.
# X: n by p matrix of the design matrix.
# L: p by p matrix of the penalty weight matrix.
# lambda.L: tuning parameter of the penalty weight matrix.
# lambda.1: tuning parameter of the L_1 penalty.
# lambda.2: tuning parameter of the ridge penalty.
# normalize.L: binary variable indicating whether the penalty weight matrix 
# needs to be normalized beforehand.
# K: number of folds in cross-validation.
# sigma.error: error standard deviation. If NULL, scale lasso is applied.
# enable.group.test: binary variable indicating whether group tests should be enabled.
# eta, C: parameters of the grace test; see Zhao & Shojaie (2016) for reference.
# ----------------------------------------------------------------------------
# Outputs:
# intercept: intercept of the linear regression model.
# beta: regression coefficient of the linear regression model.
# pvalue: p-value of individual hypothesis tests.
# group.test: function to perform group-wise hypothesis tet.

grace.test <- function(Y, X, L = NULL, lambda.L = NULL, lambda.2 = 0, normalize.L = FALSE, eta = 0.05, C = 4 * sqrt(3), K = 10, sigma.error = NULL){
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
  
  if(is.null(sigma.error)){
    sig.L <- scalreg(X, Y)$hsigma   # Error standard deviation from the scaled lasso
  }else{
    sig.L <- sigma.error
  }
  lam <- sig.L * C * sqrt(log(p) / n) / 2 # Lasso initial tuning parameter
  beta.init <- glmnet(X, Y, lambda = lam, intercept = FALSE)$beta[, 1]   # Initial estimator
  covmatrix <- solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% t(X) %*% X %*% solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p))
  se <- sig.L * sqrt(diag(covmatrix)) # Standard error
  bias <- solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% (lambda.L * L + lambda.2 * diag(p)) %*% beta.init  # Bias correction
  targ <- abs(solve(t(X) %*% X + lambda.L * L + lambda.2 * diag(p)) %*% (lambda.L * L + lambda.2 * diag(p)))
  diag(targ) <- 0
  correct <- apply(targ, 1, max) * (log(p) / n)^(0.5 - eta)  # Gamma
  
  teststat <- betahat + bias   # Bias corrected test statistic
  Tstat <- (abs(teststat) - correct) * (abs(teststat) > correct)
  pval <- 2 * (1 - pnorm(abs(Tstat / se)))
  
  group.testing.function <- function(group){
    if(!is.vector(group)){
      stop("Error: you need to supply a vector of variable indices of the group.")
    }else if(length(group) < 2){
      stop("Error: the size of the group needs to be at least two.")
    }
    test.stat = max(abs(Tstat[group] / se[group]))
    subcovmatrix <- diag(1 / sqrt(diag(covmatrix[group, group]))) %*% covmatrix[group, group] %*% diag(1 / sqrt(diag(covmatrix[group, group])))
    absZ <- abs(mvrnorm(n = 100000, mu = rep(0, length(group)), Sigma = subcovmatrix))
    maxZ <- apply(absZ, 1, max)
    groupp <- 1 - mean(test.stat > maxZ)
    return(groupp)
  }

  return(list(intercept = truealphahat, beta = truebetahat, pvalue = pval, group.test = group.testing.function))
}
