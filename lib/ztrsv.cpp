//
//  ztrsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trsv.h"

using LATL::TRSV;

int ztrsv_(char& uplo, char& trans, char& diag, int &n, complex<double> *A, int &ldA, complex<double> *x, int& incx)
{
   int info=-TRSV<double>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("ZTRSV ",info);
   return 0;
}
