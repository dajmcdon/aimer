#include <RcppArmadillo.h>
#include <iostream>
#include "irlb.h"

// [[Rcpp::depends(RcppArmadillo)]]

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


struct OUT{
    int nCovBest;
    int nCompBest;
    arma::colvec ncomps;
    arma::colvec ncovs;
    arma::mat cvMSE;
    arma::colvec mse;
};

struct TrainTest{
    arma::mat Xtrain;
    arma::mat Xtest;
    arma::colvec ytrain;
    arma::colvec ytest;
};



//returns a TrainTest (defined above) where the test data is the rows
//from start to stop including stop, and the train data is all of the rest.
TrainTest ttsplit(arma::mat X, arma::colvec y, int start, int stop){
    TrainTest output;
    output.Xtest = X.rows(start, stop);
    output.ytest = y.subvec(start, stop);
    if(start == 0){
        output.Xtrain = X.rows(stop + 1, X.n_rows - 1);
        output.ytrain = y.subvec(stop + 1, y.n_elem - 1);
    }
    else if(stop == y.n_elem - 1){
        output.Xtrain = X.rows(0, start - 1);
        output.ytrain = y.subvec(0, start - 1);
    }
    else{
        output.Xtrain = arma::join_cols(X.rows(0, start - 1),
                                        X.rows(stop + 1, X.n_rows - 1));
        output.ytrain = arma::join_cols(y.subvec(0, start-1), y.subvec(stop + 1, X.n_rows - 1));
    }
    return output;
}

//not the same as R's quantile function
double findThresh(arma::colvec t, int num){
    t = arma::sort(t);
    int index = t.n_elem - num;
    if(index == 0){
        return t[index] - 1; //everything will be above the threshold
    }
    return t[index - 1]; //num elements will be above the threshold
}

// [[Rcpp::export]]
Rcpp::List findThresholdAIMER(arma::mat X, arma::colvec y, arma::colvec ncomps,
                       arma::colvec nCovs,
                       int nthresh,
                       int kfold){
    arma::cube CVmse = arma::cube(nthresh, ncomps.n_elem, kfold);
    int start = 0;
    int stop;
    double foldSize = ((double) X.n_rows)/((double) kfold);
    for(int k = 0; k < kfold; k++){   //for each fold
        stop = (k + 1)*foldSize - 1;
        TrainTest tt = ttsplit(X, y, start, stop);
        arma::colvec tStats = arma::abs(marginalRegressionTT(tt.Xtrain, tt.ytrain));
        for(int i = 0; i < nthresh; i++){   //tries different number of nCovs[i] highest marginal correlations to accept
            double threshold = findThresh(tStats, nCovs[i]);                 //TODO: (after trying partial svd to speed up) move declarations outside of for loop
            MatrixInteger parti = partition(tt.Xtrain, tStats, threshold);  //      -add comments for for loops
            arma::mat Xnew = parti.matrix;
            arma::mat F = arma::trans(Xnew) * Xnew.cols(0, parti.num - 1);
            if(F.n_cols > F.n_rows){
                return Rcpp::List::create(Rcpp::Named("Error") = "F has more columns than rows",
                                          Rcpp::Named("columns") = F.n_cols,
                                          Rcpp::Named("rows") = F.n_rows);
            }
            for(int j = 0; j < ncomps.n_elem; j++){    //tries different number of ncomps[j] largest singular values to use
                arma::mat UF = arma::zeros<arma::mat>(F.n_rows, ncomps[j] + 7);
                arma::vec SF = arma::zeros<arma::vec>(ncomps[j] + 7);
                arma::mat VF = arma::zeros<arma::mat>(ncomps[j] + 7, F.n_cols);
                VF.col(0) = arma::randn(VF.n_rows);
                if((F.n_cols < (ncomps[j] + 7)) || (ncomps[j] > (0.5*F.n_cols))){
                    arma::svd(UF, SF, VF, F);                        //full svd
                }
                else{
                    int isError = irlb(F.memptr(), F.n_rows, F.n_cols, ncomps[j], SF.memptr(), UF.memptr(), VF.memptr());   //partial svd
                    if(isError != 0){
                        return Rcpp::List::create(Rcpp::Named("Error") = isError,
                                                  Rcpp::Named("ncomps") = ncomps[j],
                                                  Rcpp::Named("Fcolumns") = F.n_cols);
                    }
                }
                arma::mat Vd = UF.cols(0, ncomps[j] - 1);
                arma::vec Sd;
                if(ncomps[j] < SF.n_elem){
                    Sd = arma::sqrt(SF.subvec(0, ncomps[j] - 1));
                }
                else{
                    Sd = arma::sqrt(SF);
                }
                arma::mat VSInv = arma::zeros<arma::mat>(Vd.n_rows, Sd.n_elem);
                double sinv;
                for(int l = 0; l < VSInv.n_cols; l++){    //efficiently multiplys a diagonal matrix
                    sinv = 1/Sd(l);
                    for(int m = 0; m < VSInv.n_rows; m++){
                        VSInv(m, l) = sinv * Vd(m, l);
                    }
                }
                arma::mat Ud = Xnew * VSInv;
                arma::colvec beta = VSInv * trans(Ud) * tt.ytrain;
                arma::mat testX = partition(tt.Xtest, tStats, threshold).matrix;
                arma::colvec Yhat = testX * beta;
                CVmse(i, j, k) = arma::mean(arma::square((tt.ytest - Yhat)));
            }
        }
        start = stop + 1;
    }
    arma::mat mse = arma::mat(nthresh, ncomps.n_elem);
    int bestnthresh = 0;
    int bestncomp = 0;
    arma::vec temp;
    temp = CVmse.tube(0,0);
    double minMSE = arma::mean(temp);
    for(int i = 0; i < nthresh; i++){     //finds the minimal mse
        for(int j = 0; j < ncomps.n_elem; j++){
            temp = CVmse.tube(i, j);
            mse(i, j) = arma::mean(temp);
            if(mse(i,j) < minMSE){
                bestnthresh = i;
                bestncomp = j;
                minMSE = mse(i, j);
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("nCov.best") = nCovs[bestnthresh],
                              Rcpp::Named("ncomp.best") = ncomps[bestncomp],
                              Rcpp::Named("ncomps") = ncomps,
                              Rcpp::Named("nCovs") = nCovs,
                              Rcpp::Named("CVmse") = CVmse,
                              Rcpp::Named("mse") = mse);
}



//[[Rcpp::export]]
Rcpp::List findThresholdSel(arma::mat X, arma::colvec y, arma::colvec ncomps,
                              arma::colvec nCovs,
                              int nthresh,
                              int kfold,
                              arma::colvec nCovsSelect,
                              int nthreshSelect){
    arma::field<arma::colvec> CVmse = arma::field<arma::colvec>(nthreshSelect, ncomps.n_elem, nthresh);
    int start = 0;
    int stop;
    double foldSize = ((double) X.n_rows)/((double) kfold);
    for(int k = 0; k < kfold; k++){    //for each fold
        stop = (k + 1)*foldSize - 1;
        TrainTest tt = ttsplit(X, y, start, stop);
        arma::colvec tStats = arma::abs(marginalRegressionTT(tt.Xtrain, tt.ytrain));
        for(int i = 0; i < nthresh; i++){   //tries different number nCovs[i] highest marginal correlations to use
            double threshold = findThresh(tStats, nCovs[i]);
            MatrixInteger parti = partition(tt.Xtrain, tStats, threshold);
            arma::mat Xnew = parti.matrix;
            arma::mat F = arma::trans(Xnew) * Xnew.cols(0, parti.num - 1);
            if(F.n_cols > F.n_rows){
                return Rcpp::List::create(Rcpp::Named("Error") = "F has more columns than rows",
                                          Rcpp::Named("columns") = F.n_cols,
                                          Rcpp::Named("rows") = F.n_rows);
            }
            for(int j = 0; j < ncomps.n_elem; j++){   //tries different number ncomps[j] largest singular values to use
                arma::mat UF = arma::zeros<arma::mat>(F.n_rows, ncomps[j] + 7);
                arma::vec SF = arma::zeros<arma::vec>(ncomps[j] + 7);
                arma::mat VF = arma::zeros<arma::mat>(ncomps[j] + 7, F.n_cols);
                VF.col(0) = arma::randn(VF.n_rows);
                if((F.n_cols < (ncomps[j] + 7)) || (ncomps[j] > (0.5*F.n_cols))){
                    arma::svd(UF, SF, VF, F);                        //full svd
                }
                else{
                    int isError = irlb(F.memptr(), F.n_rows, F.n_cols, ncomps[j], SF.memptr(), UF.memptr(), VF.memptr());   //partial svd
                    if(isError != 0){
                        return Rcpp::List::create(Rcpp::Named("Error") = isError,
                                                  Rcpp::Named("ncomps") = ncomps[j],
                                                  Rcpp::Named("Fcolumns") = F.n_cols);
                    }
                }
                arma::mat Vd = UF.cols(0, ncomps[j] - 1);
                arma::vec Sd;
                if(ncomps[j] < SF.n_elem){
                    Sd = arma::sqrt(SF.subvec(0, ncomps[j] - 1));
                }
                else{
                    Sd = arma::sqrt(SF);
                }
                arma::mat VSInv = arma::zeros<arma::mat>(Vd.n_rows, Sd.n_elem);
                double sinv;
                for(int l = 0; l < VSInv.n_cols; l++){    //efficiently multiplies a diagonal matrix
                    sinv = 1/Sd(l);
                    for(int m = 0; m < VSInv.n_rows; m++){
                        VSInv(m, l) = sinv * Vd(m, l);
                    }
                }
                arma::mat Ud = Xnew * VSInv;
                arma::colvec initialBeta = VSInv * trans(Ud) * tt.ytrain;
                //Select process
                for(int l = 0; l < nthreshSelect; l++){    //tries different number nCovsSelect(l) largest betas to keep
                    arma::colvec beta = initialBeta;
                    double bthresh = findThresh(arma::abs(beta), nCovsSelect(l));
                    for(int m = 0; m < beta.n_elem; m++){   //sets the remaining betas to zero
                        if(beta(m) < bthresh && -1*beta(m) < bthresh){     //abs() only worked for integers
                            beta(m) = 0;
                        }
                    }
                    //no straightforward way to boolean index, so I use the whole test matrix and beta vector with zeros.
                    arma::mat testX = partition(tt.Xtest, tStats, threshold).matrix;
                    arma::colvec Yhat = testX * beta;
                    if(k == 0){
                        CVmse(l, j, i) = arma::colvec(kfold);
                    }
                    CVmse(l, j, i)(k) = arma::mean(arma::square((tt.ytest - Yhat)));
                }
            }
        }
        start = stop + 1;
    }
    arma::cube mse = arma::cube(nthreshSelect, ncomps.n_elem, nthresh);
    int bestnthresh = 0;
    int bestncomp = 0;
    int bestnthreshSelect = 0;
    arma::vec temp = CVmse(0, 0, 0);
    double minMSE = arma::mean(temp);
    for(int i = 0; i < nthreshSelect; i++){      //finds the smallest mse (and the best parameters)
        for(int j = 0; j < ncomps.n_elem; j++){
            for(int k = 0; k < nthresh; k++){
                temp = CVmse(i, j, k);
                mse(i, j, k) = arma::mean(temp);
                if(mse(i,j,k) < minMSE){
                    bestnthreshSelect = i;
                    bestncomp = j;
                    bestnthresh = k;
                    minMSE = mse(i, j, k);
                }
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("nCov.select.best") = nCovsSelect[bestnthreshSelect],
                              Rcpp::Named("ncomp.best") = ncomps[bestncomp],
                              Rcpp::Named("nCov.best") = nCovs[bestnthresh],
                              Rcpp::Named("ncomps") = ncomps,
                              Rcpp::Named("nCovs") = nCovs,
                              Rcpp::Named("nCovs.select") = nCovsSelect,
                              //Rcpp::Named("CVmse") = CVmse,       //Rcpp can't return this
                              Rcpp::Named("mse") = mse);
}



// [[Rcpp::export]]
arma::colvec AIMER(arma::mat X, arma::colvec y,
                   double t, double b, int d){
    arma::colvec xt = arma::abs(marginalRegressionTT(X, y));
    MatrixInteger parti = partition(X, xt, t);
    arma::mat Xnew = parti.matrix;
    arma::mat F = arma::trans(Xnew) * Xnew.cols(0, parti.num - 1);
    if(F.n_cols > F.n_rows){
        return arma::zeros<arma::colvec>(1);  //error
    }
    arma::mat UF = arma::zeros<arma::mat>(F.n_rows, d + 7);
    arma::vec SF = arma::zeros<arma::vec>(d + 7);
    arma::mat VF = arma::zeros<arma::mat>(d + 7, F.n_cols);
    VF.col(0) = arma::randn(VF.n_rows);
    if((F.n_cols < (d + 7)) || (d > (0.5*F.n_cols))){
        arma::svd(UF, SF, VF, F);                        //full svd
    }
    else{
        int isError = irlb(F.memptr(), F.n_rows, F.n_cols, d, SF.memptr(), UF.memptr(), VF.memptr());   //partial svd
        if(isError != 0){
            return arma::zeros<arma::colvec>(1); //error
        }
    }
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
    for(int i = 0; i < VSInv.n_cols; i++){   //efficiently multiplies diagonal matrix
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
