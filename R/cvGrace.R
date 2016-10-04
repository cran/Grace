# This function performs cross-validation for Grace
# Author: Sen Zhao
# Email: sen-zhao@sen-zhao.com

cvGrace <- function(X, Y, L, lambda.L, lambda.1, lambda.2, K = 10){
  lambda.1 <- unique(sort(lambda.1, decreasing = TRUE))
  lambda.L <- unique(sort(lambda.L, decreasing = TRUE))
  lambda.2 <- unique(sort(lambda.2, decreasing = TRUE))
  
  oneL1 <- (length(lambda.1) == 1)
  noL1 <- oneL1 & (sum(lambda.1 == 0) == 1)
  
  if(oneL1){
    lambda.1 <- c(lambda.1 + 0.01, lambda.1)
  }
  p <- ncol(X)
  n <- nrow(X)
  lam1 <- matrix(nrow = length(lambda.L), ncol = length(lambda.2))
  ERR <- matrix(NA, nrow = length(lambda.L), ncol = length(lambda.2))
  for(iL in 1:length(lambda.L)){
    lL <- lambda.L[iL]
    for(i2 in 1:length(lambda.2)){
      l2 <- lambda.2[i2]
      Lnew <- lL * L + l2 * diag(p)
      eL <- eigen(Lnew)
      S <- eL$vectors %*% sqrt(diag(eL$values))
      l1star <- lambda.1
      Xstar <- rbind(X, t(S)) / sqrt(2)
      Ystar <- c(Y, rep(0, p))
      gammastar <- l1star / sqrt(2) / 2 / (n + p)
      cvres <- cv.glmnet(Xstar, Ystar, lambda = gammastar, intercept = FALSE, standardize = FALSE, nfolds = K)
      lam1[iL, i2] <- which.min(abs(gammastar - cvres$lambda.min))
      if(oneL1){
        if(length(cvres$cvm) == 1){
          ERR[iL, i2] <- cvres$cvm[1]
        }else{
          ERR[iL, i2] <- cvres$cvm[2]
        }
      }else{
        ERR[iL, i2] <- cvres$cvm[lam1[iL, i2]]
      }
    }
  }
  idx <- which(ERR == min(ERR), arr.ind = TRUE)
  if(noL1){
    resl1 <- 0
  }else if(oneL1){
    resl1 <- lambda.1[2]  
  }else{
    resl1 <- lambda.1[lam1[idx]][1]
  }
  reslL <- lambda.L[idx[1, 1]]
  resl2 <- lambda.2[idx[1, 2]]
  return(c(reslL, resl1, resl2))
}