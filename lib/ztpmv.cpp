//
//  ztpmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tpmv.h"

using LATL::TPMV;

int ztpmv_(char& uplo, char& trans, char& diag, int &n, complex<double> *A, complex<double> *x, int& incx)
{
   int info=-TPMV<double>(uplo,trans,diag,n,A,x,incx);
   if(info>0)
      xerbla_("ZTPMV ",info);
   return 0;
}
