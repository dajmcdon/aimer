#' Predict Supervised PCA
#'
#' use this function to predict response for feature data.
#'
#' @param obj an object of class 'supervisedPCA' as output by approximateSuperPCA.
#' @param newx an optional matrix of points for predictions, if missing, then uses
#' original design matrix.
#' @param PCAmethod required.
#'
#' @return predictions for newx design matrix
#' @export
predictSuperPCA <- function(obj, newx = NULL,
                            PCAmethod=c("SPC", "SPC.lasso", "AIMER0", "AIMER")){
  if(is.null(newx)) newx <- obj$x  # can be used to compute fitted values
  switch (PCAmethod,
          SPC = {
            smallx <- newx[ , obj$keep, drop=FALSE]
            smallx <- scale(smallx, center=obj$xmeans, scale=FALSE)
            ## Note: at this point, their code does some screwy scalings (superpc.predict ~30).
            ## Check results later.
            yhat <- obj$ymean + smallx %*% obj$betahat
            return(yhat)
          },
          SPC.lasso = {
            smallx <- newx[ , obj$keep, drop=FALSE]
            yhat <- predict(obj$cvfit, newx=smallx, s="lambda.min")
            return(yhat)
          },
          AIMER0 = {
            x.new <- cbind(newx[ , obj$keep, drop=FALSE], newx[ , !obj$keep, drop=FALSE])
            x.new <- scale(x.new, center=obj$xmeans, scale=FALSE)
            yhat <- obj$ymean + x.new %*% obj$betahat
            return(yhat)
          },
          AIMER = {
            x.new <- cbind(newx[ , obj$keep, drop=FALSE], newx[ , !obj$keep, drop=FALSE])
            x.new <- scale(x.new, center=obj$xmeans, scale=FALSE)
            yhat <- obj$ymean + x.new[, obj$b.keep, drop=FALSE] %*% obj$betahat
            return(yhat)
          }
  )
}
