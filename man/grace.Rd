\name{grace}
\alias{grace}
\alias{grace-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Graph-Constrained Estimation and Hypothesis Testing
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculate coefficient estimates of Grace based on methods described in Li and Li (2008) and test for the association between Y and each column of X based on Zhao and Shojaie (2015).
}
\usage{
  grace(Y, X, L, lambda.L, lambda.L2 = 0, normalize.L = FALSE, 
  test = FALSE, eta = 0.05, C = 4)
}
\arguments{
  \item{Y}{
	  outcome variable.
  }
  \item{X}{
  	matrix of predictors.
  }
  \item{L}{
	  penalty weight matrix L.
  }
  \item{lambda.L}{
	  tuning parameter value for the penalty induced by the L matrix (see details).
  }
  \item{lambda.L2}{
    tuning parameter value for the ridge penalty (see details).
  }
  \item{normalize.L}{
	  whether the penalty weight matrix L should be normalized.
  }
  \item{test}{
  	whether the grace test should be performed. If not, this function will only calculate the coefficient estimates.
  }
  \item{eta}{
    sparsity parameter (see details).
  }
  \item{C}{
    parameter for the initial estimator (see details).
  }
}
\details{
The Grace estimator is defined as

\deqn{\beta_{hat} = (X^TX+lambda.L*L+lambda.L2*I)^{-1}X^TY}

In the formulation, L is the penalty weight matrix. Tuning parameters lambda.L and lambda.L2 are usually chosen by cross-validation. In practice, X and Y are standardized and centered, respectively, before estimating \eqn{\beta_hat}. The resulting estimate is then rescaled back into the original scale.

To perform the Grace test, the lasso initial estimator is calculated using lasso tuning parameter \eqn{\sigma_\epsilon C\sqrt{log(p)/n}}. The Grace test requires C to be larger than \eqn{2\sqrt{2}}. \eqn{\sigma_\epsilon} is the error standard deviation, which is estimated using the scaled lasso (Sun and Zhang 2012).

eta<0.5 is the sparsity parameter of \eqn{\beta}, which controls the level of bias correction. This parameter is usually unknown. Using larger values of eta leads to more conservative results.
}
\value{
  An R `list' with elements:
  \item{intercept}{fitted intercept.}
  \item{beta}{fitted regression coefficients.}
  \item{pvalue}{p-value based on the Grace test.}
}
\references{
  Li, C., and Li, H. (2008). Network-constrained regularization and variable selection for analysis of genomic data. \eqn{Bioinformatics}, 24, 1175-1182.
  
  Sun, T., and Zhang, C.-H. (2012). Scaled sparse linear regression. \eqn{Biometrika}, 99, 879-898.
  
  Zhao, S. and Shojaie, A. (2015). A signifiance test for graph-constrained estimation. \eqn{Biometrics}, to appear.
}
\author{
  Sen Zhao \email{senz@uw.edu}
}

\examples{

set.seed(120)
n <- 100
p <- 100

L <- matrix(0, nrow = p, ncol = p)
for(i in 1:10){
	L[((i - 1) * p / 10 + 1), ((i - 1) * p / 10 + 1):(i * (p / 10))] <- -1
}
diag(L) <- 0
ind <- lower.tri(L, diag = FALSE)
L[ind] <- t(L)[ind]
diag(L) <- -rowSums(L)

beta <- eigen(L)$vectors[, p] * 10

Sigma <- solve(L + 0.1 * diag(p))
sigma.error <- sqrt(t(beta) \%*\% Sigma \%*\% beta) / 10

X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
Y <- X \%*\% beta + rnorm(n, sd = sigma.error)

grace(Y, X, L, 100, test = TRUE)
}