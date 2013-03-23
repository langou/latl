//
//  ctbsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbsv.h"

using std::complex;
using LATL::tbsv;

int ctbsv_(char& uplo, char& trans, char& diag, int &n, int& k, complex<float> *A, int &ldA, complex<float> *x, int& incx)
{
   int info=-tbsv<float>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("CTBSV ",info);
   return 0;
}
