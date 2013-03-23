//
//  cgemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemv.h"

using LATL::gemv;
using std::complex;

int cgemv_(char& trans, int &m, int &n, complex<float>& alpha, complex<float> *A, int &ldA, complex<float> *x, int& incx, complex<float> &beta, complex<float> *y, int& incy)
{
   int info=-gemv<float>(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("CGEMV ",info);
   return 0;
}
