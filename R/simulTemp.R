#' Simulate some data from a sparse 2 factor model
#'
#' @param n number of observations
#' @param p number of predictors
#' @param p1 number of predictors which load onto the factors (requires 2*`p1` <= p)
#' @param lambdas relative weights of the factors
#' @param beta1 coefficient of y on the first factor (the second is calculated automatically to make the marginal correlation between those predictors and the response 0)
#' @param sig0 Standard deviation of the factor
#' @param sig1 Standard deviation of the noise
#'
#' @return A list containing:
#'
#' 1. `X`  --- n x p matrix of predictors
#' 2. `Y`  --- response vector of length n
#' 3. `theta`  --- regression coefficients on the predictors
#' 4. `beta`  --- regression coefficients on the factors (first is `beta1`)
#' 5. `U`  --- p x 2 matrix of factors
#' 6. `Lambda`  --- matrix of factor weights (as input)
#' 7. `SigXY`  --- marginal correlation between columns of X and Y
#' @export
#'
#' @examples
#' dat = factorModelSim1(100,1000,5,c(5,1),2,.1,.1)
#' @md
#' @importFrom stats rnorm
factorModelSim1 <- function(n, p, p1=2, lambdas=c(5,2),
                            beta1=1, sig0=.1, sig1=1){
  stopifnot(2*p1 <= p, p1 > 0, n > 0,
            length(beta1)==1, length(sig0)==1, length(sig1)==1)
  U = matrix(rep(c(1,-1),times=c(p1,3*p1)),ncol=2) # creates [[1,-1],[-1,-1]] for D
  U = rbind(U,matrix(0,p-2*p1,2)) # still sparse, the remaining rows have norm 0
  U = U/sqrt(2*p1) # normalize
  Lambda = diag(lambdas,2)
  D = diag(lambdas+sig0^2)
  W = U %*% sqrt(Lambda)
  w = W[p1+1,]
  beta2 = -w[1]*beta1 / w[2]
  beta = c(beta1,beta2)
  SigXY = W %*% beta
  theta = W %*% solve(D) %*% beta
  V = matrix(rnorm(2*n),ncol=2)
  X = V %*% sqrt(Lambda) %*% t(U) + sig0*matrix(rnorm(n*p),n,p)
  Y = V %*% beta + sig1*rnorm(n)
  ret = list(X=X,Y=Y,theta=theta,beta=beta,U=U,Lambda=Lambda,SigXY=SigXY)
  ret
}
