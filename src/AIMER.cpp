#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::colvec marginalRegressionTT(arma::mat X, arma::colvec y){
    int n = y.n_elem;
    int k = X.n_cols;
    double syy = 0;
    double sxy = 0;
    double sxx = 0;
    double numer = 0;
    double std = 0;
    arma::colvec output = arma::zeros<arma::colvec>(k);
    for(int j = 0; j < k; j++){
        sxx = 0;
        sxy = 0;
        for(int i = 0; i < n; i++){
            if(j == 0){
                syy = syy + y(i)*y(i);
            }
            sxy = sxy + X(i, j)*y(i);
            sxx = sxx + X(i, j)*X(i, j);
        }
        numer = sxy/sxx;
        std = sqrt(((syy/sxx) - numer*numer)/(n - 2));
        output(j) = numer/std;
    }
    return output;
}

struct MatrixInteger{
    arma::mat matrix;
    int num;
};

MatrixInteger partition(arma::mat X, arma::colvec xt, double t){
    int tail = 0;
    for(int lead = 0; lead < xt.n_elem; lead++){
        if(xt(lead) > t){
            X.swap_cols(lead, tail);
            tail++;
        }
    }
    MatrixInteger output;
    output.matrix = X;
    output.num = tail;
    return output;
}



// [[Rcpp::export]]
arma::colvec AIMER(arma::mat X, arma::colvec y,
                   double t, double b, int d){
    arma::colvec xt = arma::abs(marginalRegressionTT(X, y));
    MatrixInteger parti = partition(X, xt, t);
    arma::mat Xnew = parti.matrix;
    arma::mat F = arma::trans(Xnew) * Xnew.cols(0, parti.num - 1);
    arma::mat UF;
    arma::vec SF;
    arma::mat VF;
    arma::svd(UF, SF, VF, F);
    arma::mat Vd = UF.cols(0, d - 1);
    arma::vec Sd;
    if(d < SF.n_elem){
        Sd = arma::sqrt(SF.subvec(0, d - 1));
    }
    else{
        Sd = arma::sqrt(SF);
    }
    arma::mat VSInv = arma::zeros<arma::mat>(Vd.n_rows, Sd.n_elem);
    double sinv;
    for(int i = 0; i < VSInv.n_cols; i++){
        sinv = 1/Sd(i);
        for(int j = 0; j < VSInv.n_rows; j++){
            VSInv(j, i) = sinv * Vd(j, i);
        }
    }
    arma::mat Ud = Xnew * VSInv;
    arma::colvec beta = VSInv * trans(Ud) * y;
    for(int i = 0; i < beta.n_elem; i++){
        if(beta(i) < b && -1*beta(i) < b){     //abs() only worked for integers
            beta(i) = 0;
        }
    }
    return beta;
}
