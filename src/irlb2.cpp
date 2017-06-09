/*
 * irlb: Implicitly restarted Lanczos bidiagonalization partial SVD.
 * Copyright (c) 2016 by Bryan W. Lewis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

// #include <RcppArmadillo.h>

#include <R.h>
#define USE_RINTERNALS
#include <Rinternals.h>
#include <Rdefines.h>

#include "irlb.h"

#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Utils.h"
#include "R_ext/Parse.h"

#include "Matrix.h"
// #include "Matrix_stubs.c"


//using namespace Rcpp;

SEXP
  RNORM (int n)
  {
    char buf[4096];
    SEXP cmdSexp, cmdexpr, ans = R_NilValue;
    ParseStatus status;
    cmdSexp = PROTECT (allocVector (STRSXP, 1));
    snprintf (buf, 4095, "rnorm(%d)", n);
    SET_STRING_ELT (cmdSexp, 0, mkChar (buf));
    cmdexpr = PROTECT (R_ParseVector (cmdSexp, -1, &status, R_NilValue));
    if (status != PARSE_OK)
    {
      UNPROTECT (2);
      error ("invalid call");
    }
    for (int i = 0; i < length (cmdexpr); i++)
      ans = eval (VECTOR_ELT (cmdexpr, i), R_GlobalEnv);
    UNPROTECT (2);
    return ans;
  }



/* irlb: main computation function.
 * returns:
 *  0 on success,
 * -1 invalid dimensions,
 * -2 not converged
 * -3 out of memory
 * -4 starting vector near the null space of A
 *
 * all data must be allocated by caller, required sizes listed below
 */
int
irlb (double *A,                // Input data matrix (double case)
      //void * AS,                // input data matrix (sparse case)
      //int mult,                 // 0 -> use double *A, 1 -> use AS
      int m,                    // data matrix number of rows, must be > 3.
      int n,                    // data matrix number of columns, must be > 3.
      int nu,                   // dimension of solution
      //int work,                 // working dimension, must be > 3.
      //int maxit = 1000,                // maximum number of main iterations
      //int restart,              // 0->no, n>0 -> restarted algorithm of dimension n
      //double tol,               // convergence tolerance
      //double *scale,            // optional scale (NULL for no scale) size n * 2
      //double *shift,            // optional shift (NULL for no shift)
      //double *center,           // optional center (NULL for no center)
      // output values
      double *s,                // output singular values at least length nu
      double *U,                // output left singular vectors  m x work
      double *V)                // output right singular vectors n x work
      //int *ITER,                // ouput number of Lanczos iterations
      //int *MPROD,               // output number of matrix vector products
      //double eps,               // machine epsilon
      // working intermediate storage, sizes shown
      /*int lwork, double *V1,    // n x work
      double *U1,               // m x work
      double *W,                // m x work  input when restart > 0
      double *F,                // n
      double *B,                // work x work  input when restart > 0
      double *BU,               // work x work
      double *BV,               // work x work
      double *BS,               // work
      double *BW,               // lwork
      double *res,              // work
      double *T,                // lwork*/
      //double svtol,             // svtol limit
      //double *svratio)          // convtest extra storage vector of length work
{
    double eps = 2.220446e-16;
    int maxit = 1000;
    double tol = 0.000001;
    double svtol = tol;
    int work = nu + 7;
    //intermediate storage
    int lwork = 7 * work * (1 + work);
    double *svratio = (double *) Calloc (work, double);
    double *V1 = (double *) Calloc (n * work, double);
    double *U1 = (double *) Calloc (m * work, double);
    double *W = (double *) Calloc (m * work, double);
    double *F = (double *) Calloc (n, double);
    double *B = (double *) Calloc (work * work, double);
    double *BU = (double *) Calloc (work * work, double);
    double *BV = (double *) Calloc (work * work, double);
    double *BS = (double *) Calloc (work, double);
    double *BW = (double *) Calloc (lwork, double);
    double *res = (double *) Calloc (work, double);
    double *T = (double *) Calloc (lwork, double);
  double d, S, R, alpha, beta, R_F, SS;
  double *x;
  int jj, kk;
  int converged, restart = 0;     //added restart = 0
  int info, j, k = restart;
  int inc = 1;
  int mprod = 0;
  int iter = 0;
  double Smax = 0;
  SEXP FOO;

/* Check for valid input dimensions */
  if (work < 4 || n < 4 || m < 4)
    return -1;

  //if (restart == 0)
    memset (B, 0, work * work * sizeof (double));    //moved to declaration of B
/* Main iteration */
  while (iter < maxit)
    {
      j = 0;
/*  Normalize starting vector */
      if (iter == 0 && restart == 0)
        {
          d = F77_NAME (dnrm2) (&n, V, &inc);
          if (d < 2 * eps)
            return -1;
          d = 1 / d;
          F77_NAME (dscal) (&n, &d, V, &inc);
          //V = arma::normalise(V);
        }
      else
        j = k;

/* optionally apply scale */
      x = V + j * n;
      /*if (scale)                   //not using scale
        {
          x = scale + n;
          memcpy (scale + n, V + j * n, n * sizeof (double));
          for (kk = 0; kk < n; ++kk)
            x[kk] = x[kk] / scale[kk];
        }*/

      /*switch (mult)                   //not using mult
        {
        case 1:
          dsdmult ('n', m, n, AS, x, W + j * m);
          break;
        default:*/
        alpha = 1;
        beta = 0;
        F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m, x,      //looks like W is a matrix and X is uninitialized?
                            &inc, &beta, W + j * m, &inc);
        //}
      mprod++;
      R_CheckUserInterrupt ();

/* optionally apply shift in square cases m = n */
      /*if (shift)      //not using shift
        {
          jj = j * m;
          for (kk = 0; kk < m; ++kk)
            W[jj + kk] = W[jj + kk] + shift[0] * x[kk];
        }*/
/* optionally apply centering */
      /*if (center)       //not using center
        {
          jj = j * m;
          beta = F77_CALL (ddot) (&n, x, &inc, center, &inc);
          for (kk = 0; kk < m; ++kk)
            W[jj + kk] = W[jj + kk] - beta;
        }*/

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
          /*switch (mult)
            {
            case 1:
              dsdmult ('t', m, n, AS, W + j * m, F);
              break;
            default:*/
            alpha = 1.0;
            beta = 0.0;
            F77_NAME (dgemv) ("t", &m, &n, &alpha, (double *) A, &m,
                                W + j * m, &inc, &beta, F, &inc);
            //}
          mprod++;
          R_CheckUserInterrupt ();
/* optionally apply shift and scale */
          /*if (shift)                             //not using shift or scale
            {
              for (kk = 0; kk < m; ++kk)
                F[kk] = F[kk] + shift[0] * W[j * m + kk];
            }
          if (scale)
            {
              for (kk = 0; kk < n; ++kk)
                F[kk] = F[kk] / scale[kk];
            }*/
          SS = -S;
          F77_NAME (daxpy) (&n, &SS, V + j * n, &inc, F, &inc);
          orthog (V, F, T, n, j + 1, 1);

          if (j + 1 < work)
            {
              R_F = F77_NAME (dnrm2) (&n, F, &inc);
              R = 1.0 / R_F;
              if (R_F < 2 * eps)        // near invariant subspace
                {
                  FOO = RNORM(n);
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
/* optionally apply scale */
              x = V + (j + 1) * n;
              /*if (scale)                      //not using scale
                {
                  x = scale + n;
                  memcpy (x, V + (j + 1) * n, n * sizeof (double));
                  for (kk = 0; kk < n; ++kk)
                    x[kk] = x[kk] / scale[kk];
                }*/
              /*switch (mult)                     //not using mult
                {
                case 1:
                  dsdmult ('n', m, n, AS, x, W + (j + 1) * m);
                  break;
                default:*/
                alpha = 1.0;
                beta = 0.0;
                F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m,
                                    x, &inc, &beta, W + (j + 1) * m, &inc);
                //}
              mprod++;
              R_CheckUserInterrupt ();
/* optionally apply shift */
              /*if (shift)                                    //not using shift
                {
                  jj = j + 1;
                  for (kk = 0; kk < m; ++kk)
                    W[jj * m + kk] = W[jj * m + kk] + shift[0] * x[kk];
                }*/
/* optionally apply centering */
             /* if (center)                                  //not using center
                {
                  jj = (j + 1) * m;
                  beta = F77_CALL (ddot) (&n, x, &inc, center, &inc);
                  for (kk = 0; kk < m; ++kk)
                    W[jj + kk] = W[jj + kk] - beta;
                }*/
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
                  FOO = RNORM(m);
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

  //*ITER = iter;
  //*MPROD = mprod;
  return (converged == 1 ? 0 : -2);
}


//cholmod_common chol_c;
/* Need our own CHOLMOD error handler */
/*void attribute_hidden
irlba_R_cholmod_error (int status, const char *file, int line,
                       const char *message)
{
  if (status < 0)
    error ("Cholmod error '%s' at file:%s, line %d", message, file, line);
  else
    warning ("Cholmod warning '%s' at file:%s, line %d", message, file, line);
}

static const R_CallMethodDef CallEntries[] = {
  {"IRLB", (DL_FUNC) & IRLB, 16},
  {NULL, NULL, 0}
};

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void
R_init_irlba (DllInfo * dll)
{
  R_registerRoutines (dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols (dll, 0);
  M_R_cholmod_start (&chol_c);
  chol_c.final_ll = 1;  */        /* LL' form of simplicial factorization */
  /* need own error handler, that resets  final_ll (after *_defaults()) : */
  /*chol_c.error_handler = irlba_R_cholmod_error;
}

void
R_unload_irlba (DllInfo * dll)
{
  M_cholmod_finish (&chol_c);
}


void
dsdmult (char transpose, int m, int n, void * a, double *b, double *c)
{
  DL_FUNC sdmult = R_GetCCallable ("Matrix", "cholmod_sdmult");
  int t = transpose == 't' ? 1 : 0;
  CHM_SP cha = (CHM_SP) a;

  cholmod_dense chb;
  chb.nrow = transpose == 't' ? m : n;
  chb.d = chb.nrow;
  chb.ncol = 1;
  chb.nzmax = chb.nrow;
  chb.xtype = cha->xtype;
  chb.dtype = 0;
  chb.x = (void *) b;
  chb.z = (void *) NULL;

  cholmod_dense chc;
  chc.nrow = transpose == 't' ? n : m;
  chc.d = chc.nrow;
  chc.ncol = 1;
  chc.nzmax = chc.nrow;
  chc.xtype = cha->xtype;
  chc.dtype = 0;
  chc.x = (void *) c;
  chc.z = (void *) NULL;

  double one[] = { 1, 0 }, zero[] = { 0, 0};
  sdmult (cha, t, one, zero, &chb, &chc, &chol_c);
}*/
