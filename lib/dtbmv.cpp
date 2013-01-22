//
//  dtbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbmv.h"

using latl::tbmv;

int dtbmv_(char& uplo, char& trans, char& diag, int &n, int& k, double *A, int &ldA, double *x, int& incx)
{
   int info=-tbmv<double>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("DTBMV ",info);
   return 0;
}
