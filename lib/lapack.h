//
//  lapack.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _lapack_h
#define _lapack_h

#include <complex>
using std::complex;

extern "C"
{
   // error handler
   
   int xerbla_(const char *name, int &info);
   
   // machine constants and NaN detection
   
   float slamch_(char &opt);
   double dlamch_(char &opt);
   int disnan_(double &din);
   int sisnan_(float &din);
   
   // scalar division routines
   
   double dlapy2_(double &x, double &y);
   float slapy2_(float &x, float &y);
   double dlapy3_(double &x, double &y, double &z);
   float slapy3_(float &x, float &y, float &z);
   int sladiv_(float &A, float &B, float &C, float &D, float &P, float &Q);
   int dladiv_(double &A, double &B, double &C, double &D, double &P, double &Q);
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
   complex<double> zladiv_(complex<double> &x, complex<double> &y);
   complex<float> cladiv_(complex<float> &x, complex<float> &y);
#else
   int zladiv_(complex<double>& ret, complex<double> &x, complex<double> &y);
   int cladiv_(complex<float>& ret, complex<float> &x, complex<float> &y);
#endif
   
   // matrix norms
   
   double dlangb_(char &norm,int &n,int &kl,int &ku,double *A,int &ldA,double *work);
   double dlange_(char &norm,int &m,int &n,double *A,int &ldA,double *work);
   double dlangt_(char &norm,int &n,double *dl,double *d,double *du);
   double dlanhs_(char &norm,int &n,double *A,int &ldA,double *work);
   double dlansb_(char &norm,char &uplo,int &n,int& k,double *A,int &ldA,double *work);
   double dlansf_(char &norm,char &trans,char &uplo,int &n,double *A,double *work);
   double dlansp_(char &norm,char &uplo,int &n,double *A,double *work);
   double dlanst_(char &norm,int &n,double *d,double *e);
   double dlansy_(char &norm,char &uplo,int &n,double *A,int &ldA,double *work);
   double dlantb_(char &norm,char &uplo,char *diag,int &n,int& k,double *A,int &ldA,double *work);
   double dlantp_(char &norm,char &uplo,char *diag,int &n,double *A,double *work);
   double dlantr_(char &norm,char &uplo,char *diag,int &m,int &n,double *A,int &ldA,double *work);
   
   float slangb_(char &norm,int &n,int &kl,int &ku,float *A,int &ldA,float *work);
   float slange_(char &norm,int &m,int &n,float *A,int &ldA,float *work);
   float slangt_(char &norm,int &n,float *dl,float *d,float *du);
   float slanhs_(char &norm,int &n,float *A,int &ldA,float *work);
   float slansb_(char &norm,char &uplo,int &n,int& k,float *A,int &ldA,float *work);
   float slansf_(char &norm,char &trans,char &uplo,int &n,float *A,float *work);
   float slansp_(char &norm,char &uplo,int &n,float *A,float *work);
   float slanst_(char &norm,int &n,float *d,float *e);
   float slansy_(char &norm,char &uplo,int &n,float *A,int &ldA,float *work);
   float slantb_(char &norm,char &uplo,char *diag,int &n,int& k,float *A,int &ldA,float *work);
   float slantp_(char &norm,char &uplo,char *diag,int &n,float *A,float *work);
   float slantr_(char &norm,char &uplo,char *diag,int &m,int &n,float *A,int &ldA,float *work);
   
   double zlangb_(char &norm,int &n,int &kl,int &ku,complex<double> *A,int &ldA,double *work);
   double zlange_(char &norm,int &m,int &n,complex<double> *A,int &ldA,double *work);
   double zlangt_(char &norm,int &n,complex<double> *dl,complex<double> *d,complex<double> *du);
   double zlanhb_(char &norm,char &uplo,int &n,int& k,complex<double> *A,int &ldA,double *work);
   double zlanhe_(char &norm,char &uplo,int &n,complex<double> *A,int &ldA,double *work);
   double zlanhf_(char &norm,char &trans,char &uplo,int &n,complex<double> *A,double *work);
   double zlanhp_(char &norm,char &uplo,int &n,complex<double> *A,double *work);
   double zlanhs_(char &norm,int &n,complex<double> *A,int &ldA,double *work);
   double zlanht_(char &norm,int &n,complex<double> *d,double *e);
   double zlansb_(char &norm,char &uplo,int &n,int& k,complex<double> *A,int &ldA,double *work);
   double zlansp_(char &norm,char &uplo,int &n,complex<double> *A,double *work);
   double zlansy_(char &norm,char &uplo,int &n,complex<double> *A,int &ldA,double *work);
   double zlantb_(char &norm,char &uplo,char *diag,int &n,int& k,complex<double> *A,int &ldA,double *work);
   double zlantp_(char &norm,char &uplo,char *diag,int &n,complex<double> *A,double *work);
   double zlantr_(char &norm,char &uplo,char *diag,int &m,int &n,complex<double> *A,int &ldA,double *work);
   
   float clangb_(char &norm,int &n,int &kl,int &ku,complex<float> *A,int &ldA,float *work);
   float clange_(char &norm,int &m,int &n,complex<float> *A,int &ldA,float *work);
   float clangt_(char &norm,int &n,complex<float> *dl,complex<float> *d,complex<float> *du);
   float clanhb_(char &norm,char &uplo,int &n,int& k,complex<float> *A,int &ldA,float *work);
   float clanhe_(char &norm,char &uplo,int &n,complex<float> *A,int &ldA,float *work);
   float clanhf_(char &norm,char &trans,char &uplo,int &n,complex<float> *A,float *work);
   float clanhp_(char &norm,char &uplo,int &n,complex<float> *A,float *work);
   float clanhs_(char &norm,int &n,complex<float> *A,int &ldA,float *work);
   float clanht_(char &norm,int &n,complex<float> *d,float *e);
   float clansb_(char &norm,char &uplo,int &n,int& k,complex<float> *A,int &ldA,float *work);
   float clansp_(char &norm,char &uplo,int &n,complex<float> *A,float *work);
   float clansy_(char &norm,char &uplo,int &n,complex<float> *A,int &ldA,float *work);
   float clantb_(char &norm,char &uplo,char *diag,int &n,int& k,complex<float> *A,int &ldA,float *work);
   float clantp_(char &norm,char &uplo,char *diag,int &n,complex<float> *A,float *work);
   float clantr_(char &norm,char &uplo,char *diag,int &m,int &n,complex<float> *A,int &ldA,float *work);

   int dlarfg_(int &n, double &alpha, double *x, int &incx, double &tau);
   int slarfg_(int &n, float &alpha, float *x, int &incx, float &tau);
   int zlarfg_(int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> &tau);
   int clarfg_(int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> &tau);
   int slacpy_(char &uplo, int &m, int &n, float *A, int &ldA, float *B, int &ldB);
   int dlacpy_(char &uplo, int &m, int &n, double *A, int &ldA, double *B, int &ldB);
   int clacpy_(char &uplo, int &m, int &n, complex<float> *A, int &ldA, complex<float> *B, int &ldB);
   int zlacpy_(char &uplo, int &m, int &n, complex<double> *A, int &ldA, complex<double> *B, int &ldB);
   int dgetrf_(int &m, int &n, double *A, int &ldA, int *IPIV, int &info);
   int sgetrf_(int &m, int &n, float *A, int &ldA, int *IPIV, int &info);
   int cgetrf_(int &m, int &n, complex<float> *A, int &ldA, int *IPIV, int &info);
   int zgetrf_(int &m, int &n, complex<double> *A, int &ldA, int *IPIV, int &info);
   int csymv_(char &uplo, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zsymv_(char &uplo, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   int csyr_(char &uplo, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *A, int &ldA);
   int zsyr_(char &uplo, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *A, int &ldA);
   int sspmv_(char &uplo, int &n, complex<float> &alpha, complex<float> *A, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zspmv_(char &uplo, int &n, complex<double> &alpha, complex<double> *A, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   int cspr_(char &uplo, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *A);
   int zspr_(char &uplo, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *A);
   int slauu2_(char &uplo, int &n, float *A, int &ldA, int &info);
   int dlauu2_(char &uplo, int &n, double *A, int &ldA, int &info);
   int clauu2_(char &uplo, int &n, complex<float> *A, int &ldA, int &info);
   int zlauu2_(char &uplo, int &n, complex<double> *A, int &ldA, int &info);
   int slauum_(char &uplo, int &n, float *A, int &ldA, int &info);
   int dlauum_(char &uplo, int &n, double *A, int &ldA, int &info);
   int clauum_(char &uplo, int &n, complex<float> *A, int &ldA, int &info);
   int zlauum_(char &uplo, int &n, complex<double> *A, int &ldA, int &info);
}
#endif
