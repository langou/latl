//
//  ctrmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmv.h"

using LATL::TRMV;
using std::complex;

int ctrmv_(char& uplo, char& trans, char& diag, int &n, complex<float> *A, int &ldA, complex<float> *x, int& incx)
{
   int info=-TRMV<float>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("CTRMV ",info);
   return 0;
}
