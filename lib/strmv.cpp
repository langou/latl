//
//  strmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmv.h"

using LATL::trmv;

int strmv_(char& uplo, char& trans, char& diag, int &n, float *A, int &ldA, float *x, int& incx)
{
   int info=-trmv<float>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("STRMV ",info);
   return 0;
}
