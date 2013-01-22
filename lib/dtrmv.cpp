//
//  dtrmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmv.h"

using latl::trmv;

int dtrmv_(char& uplo, char& trans, char& diag, int &n, double *A, int &ldA, double *x, int& incx)
{
   int info=-trmv<double>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("DTRMV ",info);
   return 0;
}
