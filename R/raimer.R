#'performs Amplified, Initially Marginal, Eigenvector Regression (AIMER) for fixed t, b, and d
#'
#'@param X required, design matrix with dimension (n, p).
#'@param y required, response vector with dimension n.
#'@param t required, threshold for marginal coefficients.
#'@param b required, threshold for output coefficients.
#'@param d required, number of principal components.
#'
#'@return coefficient vector of length p.
#'
#' @export
raimer <- function(X, y, nCovs, nCovsSelect, nComps){ #FIX VARIABLES
  if(is.na(nCovs) || is.nan(nCovs) || !is.numeric(nCovs) || nCovs < 0 || nCovs %% 1 != 0){
    stop("nCovs must be an integer greater than 0")
  }
  if(is.na(nCovsSelect) || is.nan(nCovsSelect) || !is.numeric(nCovsSelect) || nCovsSelect < 0 || nCovsSelect %% 1 != 0){
    stop("nCovsSelect must be an integer greater than 0")
  }
  if(is.na(nComps) || is.nan(nComps) || !is.numeric(nComps) || nComps < 0 || nComps %% 1 != 0){
    stop("nComps must be an integer greater than 0")
  }
  if(nComps > ncol(X)){
    stop("nComps must be less than the number of columns of X")
  }
  X = as.matrix(X)
  meanx = colMeans(X)
  y = as.vector(y)
  meany = mean(y)
  y = y - meany
  X = scale(X, scale = FALSE)
  out = list(beta = as.vector(AIMER(X, y, nCovs, nCovsSelect, nComps)))
  out$fitted = X %*% out$beta
  out$residuals = y - out$fitted
  out$meany = meany
  out$meanx = meanx
  class(out) = 'aimer'
  out
}


#' Find Optimal Threshold for Amplified, Initially Marginal,
#' Eigenvector Regression (AIMER) Without Further Selection
#'
#' find the optimal number of covariates and number of components using
#' kfold cross-validation for AIMER0
#'
#' @param x required, design matrix with dimension (n,p).
#' @param y required, response vector with dimension n.
#' @param ncomps required, number of components, can be an integer or
#' a vector of integers.
#' @param nCovs optional, a vector of possible numbers of covariates.
#' @param nCovs.min optional, the smallest number of covariates, default
#' as max(\code{ncomps})+2.
#' @param nCovs.max optional, the largest number of covariates, default
#' as number of rows of \code{x}.
#' @param nthresh optional, how many \code{nCovs} to be tested, default as 25.
#' @param kfold required, the number of k in kfold cross-validation,
#' default as 10.
#' 
#' @return an object of class 'supervisedPCACV', a list with the following components
#' \item{nCov.best}{the best number of covariates which gives smallest mse
#' in cross-validation}
#' \item{ncomp.best}{the best number of components which gives smallest mse
#' in cross-validation}
#' \item{ncomps}{all the tested ncomps}
#' \item{nCovs}{all the tested numbers of covariates}
#' \item{CVmse}{mse in cross-validation with respect to each value
#' of \code{nCovs}, \code{ncomps}, \code{kfold}}
#' \item{mse}{average mse in cross-validation over k folds}
#'
#'
#' @export
findThresholdAIMER0 <- function(X, y, nComps, nCovs = NULL, 
                                nCovs.min = ifelse(is.null(nCovs), max(nComps)+2, min(nCovs)),
                                nCovs.max = ifelse(is.null(nCovs), nrow(X), max(nCovs)),
                                nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                                kfold = 10) {
    X = as.matrix(X)
    y = as.vector(y)
    y = y - mean(y)
    X = scale(X, scale = FALSE)
    indeces = sample(1:nrow(X), nrow(X))
    X = X[indeces,]
    y = y[indeces]
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    out = findThresholdAIMER(X, y, nComps, nCovs, nthresh, kfold)
    class(out) = 'supervisedPCACV'
    out$ncomps = as.vector(out$ncomps)
    out$nCovs = as.vector(out$nCovs)
    return(out)
}


#' Find Optimal Threshold for Amplified, Initially Marginal,
#' Eigenvector Regression (AIMER) With Further Selection
#'
#' find the optimal number of covariates, number of components, and number of covariates
#' in selection process for AIMER method using kfold cross-validation.
#'
#' @param x required, design matrix with dimension (n,p).
#' @param y required, response vector with dimension n.
#' @param ncomps required, number of components, can be an integer or
#' a vector of integers.
#' @param nCovs optional, a vector of possible numbers of covariates.
#' @param nCovs.min optional, the smallest number of covariates, default
#' as max(\code{ncomps})+2.
#' @param nCovs.max optional, the largest number of covariates, default
#' as number of rows of \code{x}.
#' @param nthresh optional, how many \code{nCovs} to be tested, default as 25.
#' @param nCovs.select optional, a vector of possible numbers of covariates
#' in selection process.
#' @param nCovs.min.select optional, the smallest number of covariates
#' in selection process, default as max(\code{ncomps})+2.
#' @param nCovs.max.select optional, the largest number of covariates
#' in selection process, default as number of rows of \code{x}.
#' @param nthresh.select optional, how many \code{nCovs.select} to
#' be tested in selection process, default as 25.
#' @param kfold required, the number of k in kfold cross-validation,
#' default as 10.
#'
#'
#' @return an object of class 'supervisedPCACV', a list with the following components
#' \item{nCov.select.best}{the best number of covariates in selection process
#' which gives smallest mse in cross-validation}
#' \item{nCov.best}{the best number of covariates which gives smallest mse
#' in cross-validation}
#' \item{ncomp.best}{the best number of components which gives smallest mse
#' in cross-validation}
#' \item{nCovs.select}{all the tested numbers of covariates in selection process}
#' \item{ncomps}{all the tested ncomps}
#' \item{nCovs}{all the tested numbers of covariates}
#' \item{mse}{average mse in cross-validation over k folds}
#'
#' @export
findThresholdSelect <- function (X, y, ncomps, nCovs = NULL,
                                  nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)),
                                  nCovs.max = ifelse(is.null(nCovs), nrow(X), max(nCovs)),
                                  nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                                  nCovs.select = NULL,
                                  nCovs.min.select = ifelse(is.null(nCovs.select), max(ncomps)+2, min(nCovs.select)),
                                  nCovs.max.select = ifelse(is.null(nCovs.select), nrow(X), max(nCovs.select)),
                                  nthresh.select = ifelse(is.null(nCovs.select), 25, length(nCovs.select)),
                                  kfold = 10){
    X = as.matrix(X)
    meanx = colMeans(X)
    y = as.vector(y)
    meany = mean(y)
    y = y - meany
    X = scale(X, scale = FALSE)
    indeces = sample(1:nrow(X), nrow(X))
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    if(is.null(nCovs.select)){
        nCovs.select <- round(seq(from=nCovs.min.select, to=nCovs.max.select,
                                  length.out=nthresh.select))
    }
    out = findThresholdSel(X[indeces, ], y[indeces], ncomps, nCovs, nthresh, kfold, nCovs.select, nthresh.select)
    class(out) = 'aimerCV'
    out$ncomps = as.vector(out$ncomps)
    out$nCovs = as.vector(out$nCovs)
    out$nCovsSelect = as.vector(out$nCovsSelect)
    out$beta = as.vector(out$beta)
    out$fitted = X %*% out$beta
    out$residuals = y - out$fitted
    out$meanx = meanx
    out$meany = meany
    return(out)
}