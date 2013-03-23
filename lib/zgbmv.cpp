//
//  zgbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gbmv.h"

using LATL::GBMV;
using std::complex;

int zgbmv_(char& trans, int &m, int &n, int& kL, int& kU, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int& incx, complex<double> &beta, complex<double> *y, int& incy)
{
   int info=-GBMV<double>(trans,m,n,kL,kU,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("ZGBMV ",info);
   return 0;
}
