//
//  chpr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
#include "blas.h"
#include "hpr.h"

using latl::hpr;
using std::complex;

int chpr2_(char& uplo, int &n, complex<float> &alpha, complex<float> *x, int& incx, complex<float> *y, int& incy, complex<float> *A)
{
   int info=-hpr<float>(uplo,n,alpha,x,incx,y,incy,A);
   if(info>0)
      xerbla_("CHPR2  ",info);
   return 0;
}
