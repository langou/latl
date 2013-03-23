//
//  zhbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hbmv.h"

using std::complex;
using LATL::hbmv;

int zhbmv_(char& uplo, int &n, int& k, complex<double> &alpha,complex<double> *A, int &ldA, complex<double> *x, int& incx,complex<double> &beta, complex<double> *y, int& incy)
{
   int info=-hbmv<double>(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("ZHBMV ",info);
   return 0;
}
