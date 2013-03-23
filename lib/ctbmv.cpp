//
//  ctbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbmv.h"

using LATL::tbmv;
using std::complex;

int ctbmv_(char& uplo, char& trans, char& diag, int &n, int& k, complex<float> *A, int &ldA, complex<float> *x, int& incx)
{
   int info=-tbmv<float>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("CTBMV ",info);
   return 0;
}
