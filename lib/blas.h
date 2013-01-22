//
//  blas.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _blas_h
#define _blas_h

#include <complex>

using std::complex;

typedef size_t integer;

extern "C"
{
   int xerbla_(const char *name, int& info);
   
   int srotm_(int &n, float *x, int &incx, float *y, int &incy, float *param);
   int srotmg_(float &d1, float &d2, float &x1, float &y1, float *param);
   int drotm_(int &n, double *x, int &incx, double *y, int &incy, double *param);
   int drotmg_(double &d1, double &d2, double &x1, double &y1, double *param);
   
   int saxpy_(int &n, float& alpha, float *x, int &incx, float *y, int &incy);
   int daxpy_(int &n, double& alpha , double *x, int &incx, double *y, int &incy);
   int caxpy_(int &n, complex<float>& alpha, complex<float> *x, int &incx, complex<float> *y, int &incy);
   int zaxpy_(int &n, complex<double>& alpha, complex<double> *x, int &incx, complex<double> *y, int &incy);
   
   int isamax_(int &n, float *x, int &incx);
   int idamax_(int &n, double *x, int &incx);
   int icamax_(int &n, complex<float> *x, int &incx);
   int izamax_(int &n, complex<double> *x, int &incx);
   
   float sasum_(int &n, float *x, int &incx);
   double dasum_(int &n, double *x, int &incx);
   float scasum_(int &n, complex<float> *x, int &incx);
   double dzasum_(int &n, complex<double> *x, int &incx);
   
   int scopy_(int &n, float *x, int &incx, float *y, int &incy);
   int dcopy_(int &n, double *x, int &incx, double *y, int &incy);
   int ccopy_(int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   int zcopy_(int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
   
   float snrm2_(int &n, float *x, int &incx);
   double dnrm2_(int &n, double *x, int &incx);
   float scnrm2_(int &n, complex<float> *x, int &incx);
   double dznrm2_(int &n, complex<double> *x, int &incx);
   
   int srot_(int &n, float *x, int &incx, float *y, int &incy, float& c, float& s);
   int drot_(int &n, double *x, int &incx, double *y, int &incy, double& c, double& s);
   int csrot_(int &n, complex<float> *x, int &incx, complex<float> *y, int &incy, float& c, float& s);
   int zdrot_(int &n, complex<double> *x, int &incx, complex<double> *y, int &incy, double& c, double& s);
   
   int srotg_(float& a, float& b, float& c, float& s);
   int drotg_(double& a, double& b, double& c, double& s);
   int crotg_(complex<float>& a,complex<float>& b, float& c, complex<float>& s);
   int zrotg_(complex<double>& a,complex<double>& b, double& c, complex<double>& s);
   
   int sscal_(int &n, float &a, float *x, int &incx);
   int dscal_(int &n, double &a, double *x, int &incx);
   int cscal_(int &n, complex<float> &a, complex<float> *x, int &incx);
   int csscal_(int &n, float &a, complex<float> *x, int &incx);
   int zscal_(int &n, complex<double> &a, complex<double> *x, int &incx);
   int zdscal_(int &n, double &a, complex<double> *x, int &incx);
   
   int sswap_(int &n, float *x, int &incx, float *y, int &incy);
   int dswap_(int &n, double *x, int &incx, double *y, int &incy);
   int cswap_(int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   int zswap_(int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
   
   float sdot_(int &n, float *x, int &incx, float *y, int &incy);
   double ddot_(int &n, double *x, int &incx, double *y, int &incy);
   
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
   complex<float> cdotc_(int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   complex<float> cdotu_(int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   complex<double> zdotc_(int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
   complex<double> zdotu_(int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
#else
   int cdotc_(complex<float> &ret, int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   int cdotu_(complex<float> &ret, int &n, complex<float> *x, int &incx, complex<float> *y, int &incy);
   int zdotc_(complex<double> &ret, int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
   int zdotu_(complex<double> &ret, int &n, complex<double> *x, int &incx, complex<double> *y, int &incy);
#endif
   float sdsdot_(int &n, float *B, float *x, int &incx, float *y, int &incy);
   double dsdot_(int &n, float *x, int &incx, float *y, int &incy);
   
   int sgemv_(char &trans, int &m, int &n, float &alpha, float *A, int &ldA, float *x, int &incx, float &beta, float *y, int &incy);
   int dgemv_(char &trans, int &m, int &n, double &alpha, double *A, int &ldA, double *x, int &incx, double &beta, double *y, int &incy);
   int cgemv_(char &trans, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zgemv_(char &trans, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   
   int sgbmv_(char &trans, int &m, int &n, int &KL, int &KU, float &alpha, float *A, int &ldA, float *x, int &incx, float &beta, float *y, int &incy);
   int dgbmv_(char &trans, int &m, int &n, int &KL, int &KU, double &alpha, double *A, int &ldA, double *x, int &incx, double &beta, double *y, int &incy);
   int cgbmv_(char &trans, int &m, int &n, int &KL, int &KU, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zgbmv_(char &trans, int &m, int &n, int &KL, int &KU, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   
   int ssymv_(char &uplo, int &n, float &alpha, float *A, int &ldA, float *x, int &incx, float &beta, float *y, int &incy);
   int dsymv_(char &uplo, int &n, double &alpha, double *A, int &ldA, double *x, int &incx, double &beta, double *y, int &incy);
   int chemv_(char &uplo, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zhemv_(char &uplo, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   
   int ssbmv_(char &uplo, int &n, int &K, float &alpha, float *A, int &ldA, float *x, int &incx, float &beta, float *y, int &incy);
   int dsbmv_(char &uplo, int &n, int &K, double &alpha, double *A, int &ldA, double *x, int &incx, double &beta, double *y, int &incy);
   int chbmv_(char &uplo, int &n, int &K, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zhbmv_(char &uplo, int &n, int &K, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   
   int sspmv_(char &uplo, int &n, float &alpha, float *A, float *x, int &incx, float &beta, float *y, int &incy);
   int dspmv_(char &uplo, int &n, double &alpha, double *A, double *x, int &incx, double &beta, double *y, int &incy);
   int chpmv_(char &uplo, int &n, complex<float> &alpha, complex<float> *A, complex<float> *x, int &incx, complex<float> &beta, complex<float> *y, int &incy);
   int zhpmv_(char &uplo, int &n, complex<double> &alpha, complex<double> *A, complex<double> *x, int &incx, complex<double> &beta, complex<double> *y, int &incy);
   
   int strmv_(char &uplo, char &trans, char &diag, int &n, float *A, int &ldA, float *x, int &incx);
   int dtrmv_(char &uplo, char &trans, char &diag, int &n, double *A, int &ldA, double *x, int &incx);
   int ctrmv_(char &uplo, char &trans, char &diag, int &n, complex<float> *A, int &ldA, complex<float> *x, int &incx);
   int ztrmv_(char &uplo, char &trans, char &diag, int &n, complex<double> *A, int &ldA, complex<double> *x, int &incx);
   
   int stbmv_(char &uplo, char &trans, char &diag, int &n, int &K, float *A, int &ldA, float *x, int &incx);
   int dtbmv_(char &uplo, char &trans, char &diag, int &n, int &K, double *A, int &ldA, double *x, int &incx);
   int ctbmv_(char &uplo, char &trans, char &diag, int &n, int &K, complex<float> *A, int &ldA, complex<float> *x, int &incx);
   int ztbmv_(char &uplo, char &trans, char &diag, int &n, int &K, complex<double> *A, int &ldA, complex<double> *x, int &incx);
   
   int stpmv_(char &uplo, char &trans, char &diag, int &n, float *A, float *x, int &incx);
   int dtpmv_(char &uplo, char &trans, char &diag, int &n, double *A, double *x, int &incx);
   int ctpmv_(char &uplo, char &trans, char &diag, int &n, complex<float> *A, complex<float> *x, int &incx);
   int ztpmv_(char &uplo, char &trans, char &diag, int &n, complex<double> *A, complex<double> *x, int &incx);
   
   int strsv_(char &uplo, char &trans, char &diag, int &n, float *A, int &ldA, float *x, int &incx);
   int dtrsv_(char &uplo, char &trans, char &diag, int &n, double *A, int &ldA, double *x, int &incx);
   int ctrsv_(char &uplo, char &trans, char &diag, int &n, complex<float>  *A, int &ldA, complex<float> *x, int &incx);
   int ztrsv_(char &uplo, char &trans, char &diag, int &n, complex<double>  *A, int &ldA, complex<double> *x, int &incx);
   
   int stbsv_(char &uplo, char &trans, char &diag, int &n, int &K, float *A, int &ldA, float *x, int &incx);
   int dtbsv_(char &uplo, char &trans, char &diag, int &n, int &K, double *A, int &ldA, double *x, int &incx);
   int ctbsv_(char &uplo, char &trans, char &diag, int &n, int &K, complex<float> *A, int &ldA, complex<float> *x, int &incx);
   int ztbsv_(char &uplo, char &trans, char &diag, int &n, int &K, complex<double> *A, int &ldA, complex<double> *x, int &incx);
   
   int stpsv_(char &uplo, char &trans, char &diag, int &n, float *A, float *x, int &incx);
   int dtpsv_(char &uplo, char &trans, char &diag, int &n, double *A, double *x, int &incx);
   int ctpsv_(char &uplo, char &trans, char &diag, int &n, complex<float> *A, complex<float> *x, int &incx);
   int ztpsv_(char &uplo, char &trans, char &diag, int &n, complex<double> *A, complex<double> *x, int &incx);
   
   int sger_(int &m, int &n, float &alpha, float *x, int &incx, float *y, int &incy, float *A, int &ldA);
   int dger_(int &m, int &n, double &alpha, double *x, int &incx, double *y, int &incy, double *A, int &ldA);
   int cgerc_(int &m, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *y, int &incy, complex<float> *A, int &ldA);
   int zgerc_(int &m, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *y, int &incy, complex<double> *A, int &ldA);
   int cgeru_(int &m, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *y, int &incy, complex<float> *A, int &ldA);
   int zgeru_(int &m, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *y, int &incy, complex<double> *A, int &ldA);
   
   int ssyr_(char &uplo, int &n, float &alpha, float *x, int &incx, float *A, int &ldA);
   int dsyr_(char &uplo, int &n, double &alpha, double *x, int &incx, double *A, int &ldA);
   int cher_(char &uplo, int &n, float &alpha, complex<float> *x, int &incx, complex<float> *A, int &ldA);
   int zher_(char &uplo, int &n, double &alpha, complex<double> *x, int &incx, complex<double> *A, int &ldA);
   
   int sspr_(char &uplo, int &n, float &alpha, float *x, int &incx, float *A);
   int dspr_(char &uplo, int &n, double &alpha, double *x, int &incx, double *A);
   int chpr_(char &uplo, int &n, float &alpha, complex<float> *x, int &incx, complex<float> *A);
   int zhpr_(char &uplo, int &n, double &alpha, complex<double> *x, int &incx, complex<double> *A);
   
   int ssyr2_(char &uplo, int &n, float &alpha, float *x, int &incx, float *y, int &incy, float *A, int &ldA);
   int dsyr2_(char &uplo, int &n, double &alpha, double *x, int &incx, double *y, int &incy, double *A, int &ldA);
   int cher2_(char &uplo, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *y, int &incy, complex<float> *A, int &ldA);
   int zher2_(char &uplo, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *y, int &incy, complex<double> *A, int &ldA);
   
   int sspr2_(char &uplo, int &n, float &alpha, float *x, int &incx, float *y, int &incy, float *A);
   int dspr2_(char &uplo, int &n, double &alpha, double *x, int &incx, double *y, int &incy, double *A);
   int chpr2_(char &uplo, int &n, complex<float> &alpha, complex<float> *x, int &incx, complex<float> *y, int &incy, complex<float> *A);
   int zhpr2_(char &uplo, int &n, complex<double> &alpha, complex<double> *x, int &incx, complex<double> *y, int &incy, complex<double> *A);
   
   int sgemm_(char &transA, char &transB, int &m, int &n, int &K, float &alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC);
   int dgemm_(char &transA, char &transB, int &m, int &n, int &K, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC);
   int cgemm_(char &transA, char &transB, int &m, int &n, int &K, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC);
   int chemm_(char &side, char &uplo, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC);
   int zgemm_(char &transA, char &transB, int &m, int &n, int &K, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC);
   int zhemm_(char &side, char &uplo, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC);
   
   int ssymm_(char &side, char &uplo, int &m, int &n, float &alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC);
   int dsymm_(char &side, char &uplo, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC);
   int csymm_(char &side, char &uplo, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC);
   int zsymm_(char &side, char &uplo, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC);
   
   int ssyrk_(char &uplo, char &trans, int &n, int &K, float &alpha, float *A, int &ldA, float &beta, float *C, int &ldC);
   int dsyrk_(char &uplo, char &trans, int &n, int &K, double &alpha, double *A, int &ldA, double &beta, double *C, int &ldC);
   int csyrk_(char &uplo, char &trans, int &n, int &K, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> &beta, complex<float> *C, int &ldC);
   int cherk_(char &uplo, char &trans, int &n, int &K, float &alpha, complex<float> *A, int &ldA, float &beta, complex<float> *C, int &ldC);
   int zsyrk_(char &uplo, char &trans, int &n, int &K, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> &beta, complex<double> *C, int &ldC);
   int zherk_(char &uplo, char &trans, int &n, int &K, double &alpha, complex<double> *A, int &ldA, double &beta, complex<double> *C, int &ldC);
   
   int ssyr2k_(char &uplo, char &trans, int &n, int &K, float &alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC);
   int dsyr2k_(char &uplo, char &trans, int &n, int &K, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC);
   int csyr2k_(char &uplo, char &trans, int &n, int &K, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC);
   int cher2k_(char &uplo, char &trans, int &n, int &K, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, float &beta, complex<float> *C, int &ldC);
   int zsyr2k_(char &uplo, char &trans, int &n, int &K, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC);
   int zher2k_(char &uplo, char &trans, int &n, int &K, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, double &beta, complex<double> *C, int &ldC);
   
   int strmm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, float &alpha, float *A, int &ldA, float *B, int &ldB);
   int dtrmm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB);
   int ctrmm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB);
   int ztrmm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB);
   
   int strsm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, float &alpha, float *A, int &ldA, float *B, int &ldB);
   int dtrsm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB);
   int ctrsm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB);
   int ztrsm_(char &side, char &uplo, char &trans, char &diag, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB);
   
}
#endif
