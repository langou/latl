//
//  sgemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemv.h"

using LATL::GEMV;

int sgemv_(char& trans, int &m, int &n, float& alpha, float *A, int &ldA, float *x, int& incx, float &beta, float *y, int& incy)
{
   int info=-GEMV<float>(trans,m,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("SGEMV ",info);
   return 0;
}
