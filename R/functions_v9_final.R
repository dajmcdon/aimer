#' Estimate Supervised PCA for High-dimensional Regression
#'
#' the main function, automatically chooses best number of covariates,
#' best number of component for SPC/SPC.lasso/AIMER0, and best number of
#' covariates in selection step for AIMER by minimizing the CV
#' error, returns the estimated supervised PCA object with the
#' optimal values for tuning parameters.
#'
#' @param x required, design matrix with dimension (n,p).
#' @param y required, response vector with dimension n.
#' @param ncomps required, denotes number of components, can be an integer
#' or a vector of integers.
#' @param nCovs optional, a vector of possible numbers of covariates.
#' @param nCovs.min optional, the smallest number of covariates, default as
#' max(\code{ncomps})+2.
#' @param nCovs.max optional, the largest number of covariates, default as
#' number of rows of \code{x}.
#' @param nthresh optional, how many \code{nCovs} to be tested, default as 25.
#' @param nCovs.select optional, a vector of possible numbers of covariates
#' in selection process
#' @param nCovs.min.select optional, the smallest number of covariates
#' in selection process, default as max(\code{ncomps})+2.
#' @param nCovs.max.select optional, the largest number of covariates
#' in selection process, default as number of rows of \code{x}.
#' @param nthresh.select optional, how many \code{nCovs.select} to be
#' tested in selection process, default as 25.
#' @param PCAmethod one of c("SPC", "SPC.lasso", "AIMER0", "AIMER"), defaults to SPC
#' @param kfold required, the number of k in kfold cross-validation,
#' default as 10.
#' @param progress logical, print how many folds have been calculated
#' and the the current number of covariates.
#'
#' @return two lists, the first is 'unfixTuning' of class 'supervisedPCACV'
#' which contains tuning parameter information, and the second is
#' 'fixTuning' of class 'supervisedPCA' which is the estimated model.
#' @export
estimateSuperPCA <- function(x, y, ncomps,
             nCovs = NULL,
             nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)),
             nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
             nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
             nCovs.select = NULL,
             nCovs.min.select = ifelse(is.null(nCovs.select), max(ncomps)+2, min(nCovs.select)),
             nCovs.max.select = ifelse(is.null(nCovs.select), nrow(x), max(nCovs.select)),
             nthresh.select = ifelse(is.null(nCovs.select), 25, length(nCovs.select)),
             PCAmethod = c("SPC", "SPC.lasso", "AIMER0", "AIMER"),
             kfold = 10,
             progress = FALSE){
    if(PCAmethod == "SPC" | PCAmethod == "SPC.lasso"){
        unfixTuning <- findThreshold.SPC(x=x, y=y, ncomps=ncomps, nCovs=nCovs,
                                         nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                         nthresh=nthresh,
                                         PCAmethod=PCAmethod,
                                         kfold=kfold, progress=progress)
        fixTuning <- approximateSuperPCA(x=x, y=y, ncomp=unfixTuning$ncomp.best,
                                   nCov=unfixTuning$nCov.best, PCAmethod=PCAmethod)
    } else if (PCAmethod == "AIMER0") {
        unfixTuning <- findThreshold.AIMER0(x=x, y=y, ncomps=ncomps, nCovs=nCovs,
                                        nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                        nthresh=nthresh,
                                        kfold=kfold, progress=progress)
        fixTuning <- approximateSuperPCA(x=x, y=y, ncomp=unfixTuning$ncomp.best,
                                   nCov=unfixTuning$nCov.best, PCAmethod=PCAmethod)
    } else if (PCAmethod == "AIMER") {
        unfixTuning <- findThreshold.select(x=x, y=y, ncomps=ncomps, nCovs=nCovs,
                                            nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                            nthresh=nthresh,
                                            nCovs.select=nCovs.select,
                                            nCovs.min.select=nCovs.min.select,
                                            nCovs.max.select=nCovs.max.select,
                                            nthresh.select=nthresh.select,
                                            kfold=kfold, progress=progress)
        fixTuning <- approximateSuperPCA(x=x, y=y, ncomp=unfixTuning$ncomp.best,
                                   nCov=unfixTuning$nCov.best,
                                   nCov.select=unfixTuning$nCov.select.best,
                                   PCAmethod=PCAmethod)
    }
    out <- append(unfixTuning, fixTuning)
    class(out) <- c('supervisedPCACV', 'supervisedPCA')
    return(out)
}





#' Approximate Supervised PCA for High-dimensional Regression
#'
#' performs SPC/SPC.lasso/AIMER0/AIMER for fixed \code{ncomp}, \code{nCov}
#' and \code{nCov.select}.
#'
#' @param x required, design matrix with dimension (n,p).
#' @param y required, response vector with dimension n.
#' @param ncomp required, number of components.
#' @param nCov required, number of covariates to compute SPC/SPC.lasso/AIMER0.
#' @param nCov.select required when PCAmethod="AIMER", can be < or = or > than nCov.
#' @param PCAmethod required, the method for approximate PCA.
#' @param tStats optional.
#'
#' @return a list with many components, includes much of the input, of class 'supervisedPCA'
#' \item{keep}{logical, which columns of \code{x} are used to compute SPC/AIMER0}
#' \item{tStats}{t-statistics via marginal regression for each covariate}
#' \item{gamhat}{regression coefficients on principal components}
#' \item{betahat}{regression coefficients on the kept \code{x} columns (note that order
#' has been changed)}
#' \item{ymean}{the mean of response \code{y}}
#' \item{xmeans}{the column means of design matrix \code{x}, for out-of-sample prediction}
#' \item{u}{approximated left singularvectors}
#' \item{d}{approximated sigularvalues}
#' \item{v}{approximated right sigularvectors}
#'
#' @export
approximateSuperPCA <- function (x, y, ncomp, nCov, nCov.select=NULL,
                      PCAmethod=c("SPC", "SPC.lasso", "AIMER0", "AIMER"),
                      tStats=NULL){
    if(is.null(tStats)){
        tStats <- marginalRegressionT(x=x, y=y)
    }
    threshold <- quantile(abs(tStats), 1-nCov/ncol(x))
    keep <- abs(tStats) > threshold
    switch(PCAmethod,
           SPC = {
               smallx <- x[ , keep, drop=FALSE]
               smallx.svd <- MySvd(x=smallx, ncomponent=ncomp)
               u <- smallx.svd$u
               d <- smallx.svd$d
               v <- smallx.svd$v
               xmeans <- smallx.svd$xmeans
               ymean <- mean(y)
               ## the model is y-ymean = ud %*% gam
               ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
               gamhat <- crossprod(scale(u, center=FALSE, scale=d), y-ymean) #dim(ncomp x 1)
               ## The betas are v %*% gamhat
               betahat <- v %*% gamhat #dim(keep x 1)
               out <- list(keep = keep, tStats = tStats, gamhat = gamhat, betahat = betahat,
                           ymean = ymean, xmeans = xmeans, u = u, v = v, d = d, x = x, y = y,
                           nCov = nCov, ncomp = ncomp)
           },
           SPC.lasso = {
               supervisedPCA <- approximateSuperPCA(x=x, y=y, ncomp=ncomp, nCov=nCov,
                                              PCAmethod="SPC", tStats=tStats)
               yhat <- predict(obj=supervisedPCA, newx = NULL, PCAmethod="SPC")
               cvfit <- glmnet::cv.glmnet(x=x[ , keep, drop=FALSE], y=yhat, alpha=1,
                                  intercept=TRUE, lambda.min.ratio=0.0001)
               betahat <- coef(cvfit, s = "lambda.min")
               # betahat includes intercept
               out <- list(keep = keep, tStats = tStats, gamhat = supervisedPCA$gamhat,
                           betahat = betahat, ymean = supervisedPCA$ymean,
                           xmeans = supervisedPCA$xmeans,
                           u = supervisedPCA$u, v = supervisedPCA$v, d = supervisedPCA$d,
                           x = x, y = y, nCov = nCov, ncomp = ncomp,
                           b.keep = (abs(betahat)>0), cvfit = cvfit) # b.keep includes intercept
           },
           AIMER0 = {
               x1 <- x[ , keep, drop=FALSE]
               x.new <- cbind(x1, x[ , !keep, drop=FALSE])
               L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
               svd <- MySvd(x=L.S, ncomponent=ncomp)
               v <- svd$u[ , 1:ncomp, drop=FALSE]
               d <- svd$d[1:ncomp]^(1/2) #* (p/l)^(1/4)
               u <- scale(x.new %*% v, center=FALSE, scale=d)
               print(head(x1))
               xmeans <- colMeans(x.new)
               ymean <- mean(y)
               ## the model is y-ymean = ud %*% gam
               ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
               gamhat <- crossprod(scale(u, center=FALSE, scale=d), y-ymean) #dim(ncomp x 1)
               ## The betas are v %*% gamhat
               betahat <- v %*% gamhat #dim(p x 1)
               out <- list(keep = keep, tStats = tStats, gamhat = gamhat,
                           betahat = betahat, ymean = ymean, xmeans = xmeans,
                           u = u, v = v, d = d, x = x, y = y, nCov = nCov,
                           ncomp = ncomp, b.keep = (abs(betahat)>0))
           },
           AIMER = {
               AIMER0 <- approximateSuperPCA(x=x, y=y, ncomp=ncomp, nCov=nCov,
                                   PCAmethod="AIMER0", tStats=tStats)
               b.threshold <- quantile(abs(AIMER0$betahat), 1-nCov.select/ncol(x))
               b.keep <- abs(AIMER0$betahat) > b.threshold
               out <- list(keep = keep, tStats = tStats, gamhat = AIMER0$gamhat,
                           betahat = AIMER0$betahat[b.keep],
                           ymean = AIMER0$ymean, xmeans = AIMER0$xmeans,
                           u = AIMER0$u, v = AIMER0$v, d = AIMER0$d,
                           x = x, y = y, nCov = nCov, ncomp = ncomp,
                           b.keep = b.keep, nCov.select = nCov.select)
           }
    )
    class(out) <- 'supervisedPCA'
    return(out)
}




#' Find Optimal Threshold for Supervised Principal Components
#'
#' find the optimal number of covariates and number of components using
#' kfold cross-validation for SPC and SPC.lasso
#'
#' @param x required, design matrix with dimension (n,p).
#' @param y required, response vector with dimension n.
#' @param ncomps required, denotes number of components, can be a integer or a
#' vector of integers.
#' @param nCovs optional, a vector of possible numbers of covariates.
#' @param nCovs.min optional, the smallest number of covariates, default as
#' max(\code{ncomps})+2.
#' @param nCovs.max optional, the largest number of covariates, default as
#' number of rows of \code{x}.
#' @param nthresh optional, how many \code{nCovs} to be tested, default as 25.
#' @param PCAmethod required, either SPC or SPC.lasso.
#' @param kfold required, the number of k in kfold cross-validation, default as 10
#' @param progress logical, print how many nCov/folds have been calculated
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
#' @export
findThreshold.SPC <- function (x, y, ncomps, nCovs = NULL,
                               nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)),
                               nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
                               nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                               PCAmethod = c("SPC", "SPC.lasso"),
                               kfold = 10, progress = FALSE) {
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov
    CVmse <- array(NA, dim=c(nthresh, length(ncomps), kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY)
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            smallx <- trainingX[ , keep, drop=FALSE]
            smallx.svd <- MySvd(x=smallx, ncomponent=max(ncomps))
            xmeans <- smallx.svd$xmeans
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                u <- smallx.svd$u[, 1:ncomps[j]]
                d <- smallx.svd$d[1:ncomps[j]]
                v <- smallx.svd$v[, 1:ncomps[j]]
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(keep x 1)
                if(PCAmethod == "SPC"){
                    smallx.test <- testingX[ , keep, drop=FALSE]
                    smallx.test <- scale(smallx.test, center=xmeans, scale=FALSE)
                    testingY.hat <- ymean + smallx.test %*% betahat
                    CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
                } else if (PCAmethod == "SPC.lasso") {
                    smallx.scale <- scale(smallx, center=xmeans, scale=FALSE)
                    trainingY.hat <- ymean + smallx.scale %*% betahat
                    cvfit <- glmnet::cv.glmnet(x=smallx, y=trainingY.hat, alpha=1,
                                       intercept=TRUE, lambda.min.ratio=0.0001)
                    smallx.test <- testingX[ , keep, drop=FALSE]
                    testingY.hat <- predict(cvfit, newx=smallx.test, s="lambda.min")
                    CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
                }
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')
    }
    mse <- apply(CVmse, c(1,2), mean)
    out <- list(nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[1]],
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                ncomps = ncomps, nCovs = nCovs, CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
    return(out)
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
#' @param progress logical, print how many nCov/folds have been calculated.
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
findThreshold.AIMER0 <- function (x, y, ncomps, nCovs = NULL,
                    nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)),
                    nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
                    nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                    kfold = 10, progress = FALSE) {
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov
    CVmse <- array(NA, dim=c(nthresh, length(ncomps), kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY)
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            x1 <- trainingX[ , keep, drop=FALSE]
            x.new <- cbind(x1, trainingX[ , !keep, drop=FALSE])
            L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
            svd <- MySvd(x=L.S, ncomponent=max(ncomps))
            xmeans <- colMeans(x.new)
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                v <- svd$u[ , 1:ncomps[j], drop=FALSE]
                d <- svd$d[1:ncomps[j]]^(1/2) #* (p/l)^(1/4)
                u <- scale(x.new %*% v, center=FALSE, scale=d)
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(p x 1)
                testingX.new <- cbind(testingX[ , keep, drop=FALSE], testingX[ , !keep, drop=FALSE])
                testingX.new <- scale(testingX.new, center=xmeans, scale=FALSE)
                testingY.hat <- ymean + testingX.new %*% betahat
                CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')
    }
    mse <- apply(CVmse, c(1,2), mean)
    out <- list(nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[1]],
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                ncomps = ncomps, nCovs = nCovs, CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
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
#' @param progress logical, print how many nCov/folds have been calculated.
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
#' \item{CVmse}{mse in cross-validation with respect to each value
#' of \code{nCovs}, \code{ncomps}, \code{kfold}}
#' \item{mse}{average mse in cross-validation over k folds}
#'
#' @export
findThreshold.select <- function (x, y, ncomps, nCovs = NULL,
          nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)),
          nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
          nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
          nCovs.select = NULL,
          nCovs.min.select = ifelse(is.null(nCovs.select), max(ncomps)+2, min(nCovs.select)),
          nCovs.max.select = ifelse(is.null(nCovs.select), nrow(x), max(nCovs.select)),
          nthresh.select = ifelse(is.null(nCovs.select), 25, length(nCovs.select)),
          kfold = 10, progress = FALSE) {
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    }
    if(is.null(nCovs.select)){
        nCovs.select <- round(seq(from=nCovs.min.select, to=nCovs.max.select,
                                  length.out=nthresh.select))
    }
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov, nCov.select
    CVmse <- array(NA, dim=c(nthresh.select, length(ncomps), nthresh, kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY)
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            x1 <- trainingX[ , keep, drop=FALSE]
            x.new <- cbind(x1, trainingX[ , !keep, drop=FALSE])
            L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
            svd <- MySvd(x=L.S, ncomponent=max(ncomps))
            xmeans <- colMeans(x.new)
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                v <- svd$u[ , 1:ncomps[j], drop=FALSE]
                d <- svd$d[1:ncomps[j]]^(1/2) #* (p/l)^(1/4)
                u <- scale(x.new %*% v, center=FALSE, scale=d)
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(p x 1)
                # select process
                for(l in 1:nthresh.select){
                    b.threshold <- quantile(abs(betahat), 1-nCovs.select[l]/ncol(trainingX))
                    b.keep <- abs(betahat) > b.threshold
                    testingX.new <- cbind(testingX[, keep, drop=FALSE], testingX[, !keep, drop=FALSE])
                    testingX.new <- scale(testingX.new, center=xmeans, scale=FALSE)
                    testingY.hat <- ymean + testingX.new[, b.keep, drop=FALSE] %*% betahat[b.keep]
                    CVmse[l, j, i, k] <- mean((testingY - testingY.hat)^2)
                }
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')
    }
    mse <- apply(CVmse, c(1, 2, 3), mean)
    out <- list(nCov.select.best = nCovs.select[which(mse == min(mse), arr.ind = TRUE)[1]],
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[3]],
                ncomps = ncomps, nCovs = nCovs, nCovs.select = nCovs.select,
                CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
    return(out)
}













