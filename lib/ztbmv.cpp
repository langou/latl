//
//  ztbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbmv.h"

using LATL::tbmv;
using std::complex;

int ztbmv_(char& uplo, char& trans, char& diag, int &n, int& k, complex<double> *A, int &ldA, complex<double> *x, int& incx)
{
   int info=-tbmv<double>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("ZTBMV ",info);
   return 0;
}
