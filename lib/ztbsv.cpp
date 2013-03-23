//
//  ztbsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbsv.h"

using std::complex;
using LATL::tbsv;

int ztbsv_(char& uplo, char& trans, char& diag, int &n, int& k, complex<double> *A, int &ldA, complex<double> *x, int& incx)
{
   int info=-tbsv<double>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("ZTBSV ",info);
   return 0;
}
