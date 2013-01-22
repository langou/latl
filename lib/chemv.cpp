//
//  chemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
#include "blas.h"
#include "hemv.h"

using std::complex;
using latl::hemv;

int chemv_(char& uplo, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int& incx, complex<float> &beta, complex<float> *y, int& incy)
{
   int info=-hemv<float>(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("CHEMV ",info);
   return 0;
}
