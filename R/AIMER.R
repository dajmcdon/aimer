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