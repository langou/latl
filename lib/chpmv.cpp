//
//  chpmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hpmv.h"
using LATL::hpmv;
using std::complex;

int chpmv_(char& uplo, int &n, complex<float> &alpha, complex<float> *A, complex<float> *x, int& incx, complex<float> &beta, complex<float> *y, int& incy)
{
   int info=-hpmv<float>(uplo,n,alpha,A,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("CHPMV ",info);
   return 0;
}
