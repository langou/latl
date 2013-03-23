//
//  dspmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "spmv.h"

using LATL::SPMV;

int dspmv_(char& uplo, int &n, double &alpha, double *A, double *x, int& incx, double &beta, double *y, int& incy)
{
   int info=-SPMV<double>(uplo,n,alpha,A,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("DSPMV ",info);
   return 0;
}
