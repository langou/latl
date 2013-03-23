//
//  dtrsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trsv.h"

using LATL::trsv;

int dtrsv_(char& uplo, char& trans, char& diag, int &n, double *A, int &ldA, double *x, int& incx)
{
   int info=-trsv<double>(uplo,trans,diag,n,A,ldA,x,incx);
   if(info>0)
      xerbla_("DTRSV ",info);
   return 0;
}
