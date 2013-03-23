//
//  zhpmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hpmv.h"
using LATL::hpmv;
using std::complex;

int zhpmv_(char& uplo, int &n, complex<double> &alpha, complex<double> *A, complex<double> *x, int& incx, complex<double> &beta, complex<double> *y, int& incy)
{
   int info=-hpmv<double>(uplo,n,alpha,A,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("ZHPMV ",info);
   return 0;
}
