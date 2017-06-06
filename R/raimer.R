#' @export
raimer <- function(X, y, t, b, d){
  if(is.na(t) || is.nan(t) || !is.numeric(t) || t < 0){
    stop("t must be a number greater than 0")
  }
  if(is.na(b) || is.nan(b) || !is.numeric(b) || b < 0){
    stop("b must be a number greater than 0")
  }
  if(is.na(d) || is.nan(d) || !is.numeric(d) || d < 0 || d %% 1 != 0){
    stop("d must be an integer greater than 0")
  }
  if(d > ncol(X)){
    stop("d must be less than the number of columns of X")
  }
  X = as.matrix(X)
  y = as.vector(y)
  y = y - mean(y)
  X = scale(X, scale = FALSE)
  AIMER(X, y, t, b, d)
}


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
    y = as.vector(y)
    y = y - mean(y)
    X = scale(X, scale = FALSE)
    indeces = sample(1:nrow(X), nrow(X))
    X = X[indeces,]
    y = y[indeces]
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    if(is.null(nCovs.select)){
        nCovs.select <- round(seq(from=nCovs.min.select, to=nCovs.max.select,
                                  length.out=nthresh.select))
    }
    out = findThresholdSel(X, y, ncomps, nCovs, nthresh, kfold, nCovs.select, nthresh.select)
    class(out) = 'supervisedPCACV'
    out$ncomps = as.vector(out$ncomps)
    out$nCovs = as.vector(out$nCovs)
    out$nCovsSelect = as.vector(out$nCovsSelect)
    return(out)
}