#include <RcppArmadillo.h>
#include <iostream>

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

arma::uvec randPerm(int k){
    std::srand(time(NULL));            //still doesn't quite seem random
    arma::uvec output = arma::uvec(k);
    for(int i = 0; i < k; i++){
        output[i] = i;
    }
    std::random_shuffle(output.begin(), output.end());
    return output;
}




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
double findThresh(arma::colvec t, double p){
    t = arma::sort(t);
    int index = t.n_elem*p;
    return t[index + 1];
}

// [[Rcpp::export]]
arma::uvec findThresholdAIMER(arma::mat X, arma::colvec y, arma::colvec ncomps, 
                       arma::colvec nCovs,
                       int nCovsMin,
                       int nCovsMax,
                       int nthresh,
                       int kfold, bool progress){
    arma::uvec indeces = randPerm(X.n_rows);
    X = X.rows(indeces);
    y = y.elem(indeces);
    arma::cube CVmse = arma::cube(nthresh, ncomps.n_elem, kfold);
    int start = 0;
    int stop;
    double foldSize = ((double) X.n_rows)/((double) kfold);
    for(int k = 1; k <= kfold; k++){
        stop = k*foldSize - 1;
        TrainTest tt = ttsplit(X, y, start, stop);
        arma::colvec tStats = marginalRegressionTT(tt.Xtrain, tt.ytrain);
        for(int i = 0; i < nthresh; i++){
            double threshold = findThresh(tStats, 1-(nCovs[i]/tt.Xtrain.n_cols));
            MatrixInteger parti = partition(tt.Xtrain, tStats, threshold);
            arma::mat Xnew = parti.matrix;
            arma::mat F = arma::trans(Xnew) * Xnew.cols(0, parti.num - 1);
            arma::mat UF;
            arma::vec SF;
            arma::mat VF;
            arma::svd(UF, SF, VF, F);
            for(int j = 0; j < ncomps.n_elem; j++){
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
                for(int l = 0; l < VSInv.n_cols; l++){
                    sinv = 1/Sd(l);
                    for(int m = 0; m < VSInv.n_rows; m++){
                        VSInv(m, l) = sinv * Vd(m, l);
                    }
                }
                arma::mat Ud = Xnew * VSInv;
                arma::colvec beta = VSInv * trans(Ud) * tt.ytrain;
                arma::mat testX = partition(tt.Xtest, tStats, threshold).matrix;
                arma::colvec Yhat = testX * beta;
                CVmse[i, j, k] = arma::mean(arma::square((tt.ytest - Yhat)));
            }
        }
        start = stop + 1;
    }
    return indeces;
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
