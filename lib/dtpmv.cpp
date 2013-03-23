//
//  dtpmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tpmv.h"

using LATL::tpmv;

int dtpmv_(char& uplo, char& trans, char& diag, int &n, double *A, double *x, int& incx)
{
   int info=-tpmv<double>(uplo,trans,diag,n,A,x,incx);
   if(info>0)
      xerbla_("DTPMV ",info);
   return 0;
}
