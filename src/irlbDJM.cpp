#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

void
  orthog (arma::mat& X, arma::mat& Y, arma::mat& T, int xm, int xn, int yn)
  {
    double a = 1, b = 1;
    int inc = 1;
    // T = t(X) * Y
    T = X.t() * Y;
    // Y = Y - X * T
    Y -= X * T;
  }

arma::mat  irlb (arma::mat A,                // Input data matrix (double case)
        int m,                    // data matrix number of rows, must be > 3.
        int n,                    // data matrix number of columns, must be > 3.
        int nu,                   // dimension of solution
        int work,                 // working dimension, must be > 3.
        int maxit,                // maximum number of main iterations
        int restart,              // 0->no, n>0 -> restarted algorithm of dimension n
        double tol,               // convergence tolerance
        // output values
        double *s,                // output singular values at least length nu
        double *U,                // output left singular vectors  m x work
        double *V,                // output right singular vectors n x work
        int *ITER,                // ouput number of Lanczos iterations
        int *MPROD,               // output number of matrix vector products
        double eps,               // machine epsilon
        // working intermediate storage, sizes shown
        int lwork, double *V1,    // n x work
        double *U1,               // m x work
        double *W,                // m x work  input when restart > 0
        double *F,                // n
        double *BU,               // work x work
        double *BV,               // work x work
        double *BS,               // work
        double *BW,               // lwork
        double *res,              // work
        double *T,                // lwork
        double svtol,             // svtol limit
        double *svratio)          // convtest extra storage vector of length work
  {
    double d, S, R, alpha, beta, R_F, SS;
    double *x;
    int jj, kk;
    int converged;
    int info, j, k = restart;
    int inc = 1;
    int mprod = 0;
    int iter = 0;
    double Smax = 0;
    
    /* Check for valid input dimensions */
    if (work < 4 || n < 4 || m < 4)
      return -1;
    
    
    arma::mat B(work,work);
    /* Main iteration */
    while (iter < maxit)
    {
      j = 0;
      /*  Normalize starting vector */
      if (iter == 0)
      {
        V = arma::normalise(V);
        /* d = F77_NAME (dnrm2) (&n, V, &inc);
        if (d < 2 * eps)
          return -1;
        d = 1 / d;
        F77_NAME (dscal) (&n, &d, V, &inc);
         */
      }
      else
        j = k;
      
      x = V + j * n;
      
      /* W[,j] = Ax 
      alpha = 1;
      beta = 0;
      F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m, x,
                &inc, &beta, W + j * m, &inc);
       */
      W.col(j) = 
      
      mprod++;
      R_CheckUserInterrupt ();
      
      
      if (iter > 0)
        orthog (W, W + j * m, T, m, j, 1);
      S = F77_NAME (dnrm2) (&m, W + j * m, &inc);
      if (S < 2 * eps && j == 0)
        return -4;
      SS = 1.0 / S;
      F77_NAME (dscal) (&m, &SS, W + j * m, &inc);
      
      /* The Lanczos process */
      while (j < work)
      {
        alpha = 1.0;
        beta = 0.0;
        F77_NAME (dgemv) ("t", &m, &n, &alpha, (double *) A, &m,
                  W + j * m, &inc, &beta, F, &inc);
        mprod++;
        R_CheckUserInterrupt ();
        
        SS = -S;
        F77_NAME (daxpy) (&n, &SS, V + j * n, &inc, F, &inc);
        orthog (V, F, T, n, j + 1, 1);
        
        if (j + 1 < work)
        {
          R_F = F77_NAME (dnrm2) (&n, F, &inc);
          R = 1.0 / R_F;
          if (R_F < 2 * eps)        // near invariant subspace
          {
            FOO = RNORM (n);
            for (kk = 0; kk < n; ++kk)
              F[kk] = REAL (FOO)[kk];
            orthog (V, F, T, n, j + 1, 1);
            R_F = F77_NAME (dnrm2) (&n, F, &inc);
            R = 1.0 / R_F;
            R_F = 0;
          }
          
          memmove (V + (j + 1) * n, F, n * sizeof (double));
          F77_NAME (dscal) (&n, &R, V + (j + 1) * n, &inc);
          B[j * work + j] = S;
          B[(j + 1) * work + j] = R_F;
          
          x = V + (j + 1) * n;
          
          alpha = 1.0;
          beta = 0.0;
          F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m,
                    x, &inc, &beta, W + (j + 1) * m, &inc);
          
          mprod++;
          R_CheckUserInterrupt ();
          
          /* One step of classical Gram-Schmidt */
          R = -R_F;
          F77_NAME (daxpy) (&m, &R, W + j * m, &inc, W + (j + 1) * m,
                    &inc);
          /* full re-orthogonalization of W_{j+1} */
          orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
          S = F77_NAME (dnrm2) (&m, W + (j + 1) * m, &inc);
          SS = 1.0 / S;
          if (S < 2 * eps)
          {
            FOO = RNORM (m);
            jj = (j + 1) * m;
            for (kk = 0; kk < m; ++kk)
              W[jj + kk] = REAL (FOO)[kk];
            orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
            S = F77_NAME (dnrm2) (&m, W + (j + 1) * m, &inc);
            SS = 1.0 / S;
            F77_NAME (dscal) (&m, &SS, W + (j + 1) * m, &inc);
            S = 0;
          }
          else
            F77_NAME (dscal) (&m, &SS, W + (j + 1) * m, &inc);
        }
        else
        {
          B[j * work + j] = S;
        }
        j++;
      }
      
      memmove (BU, B, work * work * sizeof (double));   // Make a working copy of B
      int *BI = (int *) T;
      F77_NAME (dgesdd) ("O", &work, &work, BU, &work, BS, BU, &work, BV,
                &work, BW, &lwork, BI, &info);
      R_F = F77_NAME (dnrm2) (&n, F, &inc);
      R = 1.0 / R_F;
      F77_NAME (dscal) (&n, &R, F, &inc);
      /* Force termination after encountering linear dependence */
      if (R_F < 2 * eps)
        R_F = 0;
      
      Smax = 0;
      for (jj = 0; jj < j; ++jj)
      {
        if (BS[jj] > Smax)
          Smax = BS[jj];
        svratio[jj] = fabs (svratio[jj] - BS[jj]) / BS[jj];
      }
      for (kk = 0; kk < j; ++kk)
        res[kk] = R_F * BU[kk * work + (j - 1)];
      /* Update k to be the number of converged singular values. */
      convtests (j, nu, tol, svtol, Smax, svratio, res, &k, &converged, S);
      if (converged == 1)
      {
        iter++;
        break;
      }
      for (jj = 0; jj < j; ++jj)
        svratio[jj] = BS[jj];
      
      alpha = 1;
      beta = 0;
      F77_NAME (dgemm) ("n", "t", &n, &k, &j, &alpha, V, &n, BV, &work, &beta,
                V1, &n);
      memmove (V, V1, n * k * sizeof (double));
      memmove (V + n * k, F, n * sizeof (double));
      
      memset (B, 0, work * work * sizeof (double));
      for (jj = 0; jj < k; ++jj)
      {
        B[jj * work + jj] = BS[jj];
        B[k * work + jj] = res[jj];
      }
      
      /*   Update the left approximate singular vectors */
      alpha = 1;
      beta = 0;
      F77_NAME (dgemm) ("n", "n", &m, &k, &j, &alpha, W, &m, BU, &work, &beta,
                U1, &m);
      memmove (W, U1, m * k * sizeof (double));
      iter++;
    }
    
    /* Results */
    memmove (s, BS, nu * sizeof (double));        /* Singular values */
    alpha = 1;
    beta = 0;
    F77_NAME (dgemm) ("n", "n", &m, &nu, &work, &alpha, W, &m, BU, &work, &beta,
              U, &m);
    F77_NAME (dgemm) ("n", "t", &n, &nu, &work, &alpha, V, &n, BV, &work, &beta,
              V1, &n);
    memmove (V, V1, n * nu * sizeof (double));
    
    *ITER = iter;
    *MPROD = mprod;
    return (converged == 1 ? 0 : -2);
  }