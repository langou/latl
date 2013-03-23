//
//  dgemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemv.h"

using LATL::gemv;

int dgemv_(char& trans, int &m, int &n, double& alpha, double *A, int &ldA, double *x, int& incx, double &beta, double *y, int& incy)
{
   int info=-gemv<double>(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("DGEMV ",info);
   return 0;
}
