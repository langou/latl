//
//  chbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hbmv.h"

using std::complex;
using latl::hbmv;

int chbmv_(char& uplo, int &n, int& k, complex<float> &alpha,complex<float> *A, int &ldA, complex<float> *x, int& incx,complex<float> &beta, complex<float> *y, int& incy)
{
   int info=-hbmv<float>(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("CHBMV ",info);
   return 0;
}
