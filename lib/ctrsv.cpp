//
//  ctrsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trsv.h"

using LATL::TRSV;

int ctrsv_(char& uplo, char& trans, char& diag, int &n, complex<float> *A, int &ldA, complex<float> *x, int& incx)
{
   int info=-TRSV<float>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("CTRSV ",info);
   return 0;
}
