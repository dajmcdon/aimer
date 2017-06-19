## Some helper functions which are not meant to be called directly

#' @export
marginalRegressionT <- function (x, y) {
  ## Purpose: gets t-statistics via marginal regression
  ## Inputs: an n x p matrix x
  ##             an n vector y
  ## Outputs: a p vector of t-statistics for the univariate (marginal regressions)
  ##                of each column of x on y
  y <- as.vector(y)
  n <- length(y)
  cx <- scale(x, scale=FALSE) # centered only
  cy <- y - mean(y)
  sxx <- colSums(cx^2)
  sxy <- crossprod(cx , cy) # t(x) %*% y but faster
  syy <- sum(cy^2)
  numer <- sxy/sxx
  std <- sqrt((syy/sxx - numer^2)/(n - 2))
  return(numer / std)
}


MySvd <- function(x, ncomponent=min(nrow(x), ncol(x))) {
  ## Purpose: performs svd quickly, if ncomponent <= min(nrow(x), ncol(x)) * .5, use irlba,
  ##          otherwise use svd directly
  ## Inputs: an n x p matrix x
  ##             an (optional) number of desired components, defaults to "all"
  ## Outputs: a list with components u, d, v, and xmeans
  ##             u is a matrix of left singular vectors of dimension n x ncomponent
  ##             v is a matrix of right singular vectors (not t(v)) of dimension p x ncomponent
  ##             d is a vector of singular values of length ncomponent
  ##             xmeans is a p-vector of the column means of x
  xmeans <- colMeans(x)
  x <- scale(x, center = xmeans, scale=FALSE)
  # I think we  should require the input of ncomp to be larger than 0

  #     if(ncomponent==0){
  #         s <- svd(x, nu=0, nv=0)
  #         return(list(d = s$d, xmeans=xmeans))
  #     }

  if (ncomponent < min(nrow(x), ncol(x)) * .5) {
    s <- irlba::irlba(A=x, nv=ncomponent, nu=ncomponent)
    u <- s$u[,1:ncomponent,drop=FALSE]
    d <- s$d[1:ncomponent]
    v <- s$v[,1:ncomponent,drop=FALSE]
  } else {
    s <- svd(x=x, nv=ncomponent, nu=ncomponent)
    u <- s$u[,1:ncomponent,drop=FALSE]
    d <- s$d[1:ncomponent]
    v <- s$v[,1:ncomponent,drop=FALSE]
  }
  return(list(u=u, d=d, v=v, xmeans=xmeans)) # v not t(v)
}



MykfoldCV <- function(x, y, k=10){
  ## Purpose: prepare training set and testing set for kfold cross-validation
  ## Inputs: an n*p matrix x as data matrix;
  ##         an n*1 vector y as response;
  ##         k is the number of folds, by default is 10.
  ## Outputs: a list of 2 components, first is the list of dataset including k lists,
  ##          each includes sets of training.x, training.y, testing.x, testing.y;
  ##          second is the folds list which indicates the index for each fold.

  # randomly select index for partition of rows of x
  n <- nrow(x)
  folds <- vector("list", k)
  breaks <- round(seq(from = 1, to = (n + 1), length = (k + 1)))
  cv.order <- sample(1:n)
  dataset <- vector("list", k)
  for(i in 1:k){
    # prepare index for the ith fold
    folds[[i]] <- cv.order[(breaks[i]):(breaks[i + 1] - 1)]
    # generate training set and testing set for ith fold
    testing.x <- x[folds[[i]], ]
    testing.y <- y[folds[[i]]]
    training.x <- x[-folds[[i]], ]
    training.y <- y[-folds[[i]]]
    dataset[[i]] <- list(testing.x = testing.x, testing.y = testing.y,
                         training.x = training.x, training.y=training.y)
  }
  out <- list(dataset=dataset, folds=folds)
  return(out)
}
