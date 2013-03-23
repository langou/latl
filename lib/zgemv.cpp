//
//  zgemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemv.h"

using LATL::gemv;
using std::complex;

int zgemv_(char& trans, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int& incx, complex<double> &beta, complex<double> *y, int& incy)
{
   int info=-gemv<double>(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("ZGEMV ",info);
   return 0;
}
