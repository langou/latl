//
//  dgbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gbmv.h"

using latl::gbmv;

int dgbmv_(char& trans, int &m, int &n, int& kL, int& kU, double &alpha, double *A, int &ldA, double *x, int& incx, double &beta, double *y, int& incy)
{
   int info=-gbmv<double>(trans,m,n,kL,kU,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("DGBMV ",info);
   return 0;
}
