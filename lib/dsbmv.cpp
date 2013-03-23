//
//  dsbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "sbmv.h"

using LATL::sbmv;

int dsbmv_(char& uplo, int &n, int& k, double &alpha, double *A, int &ldA, double *x, int& incx, double &beta, double *y, int& incy)
{
   int info=-sbmv<double>(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("DSBMV ",info);
   return 0;
}
