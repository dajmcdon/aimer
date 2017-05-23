factorModelSim1 <- function(n, p, p1, lambdas, beta1, sig0, sig1){
  U = matrix(rep(c(1,-1),times=c(p1,3*p1)),ncol=2) # creates [[1,-1],[-1,-1]] for D
  U = rbind(U,matrix(0,p-2*p1,2)) # still sparse, the remaining rows have norm 0
  U = U/sqrt(2*p1) # normalize
  Lambda = diag(lambdas)
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
