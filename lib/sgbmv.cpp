//
//  sgbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gbmv.h"

using LATL::gbmv;

int sgbmv_(char& trans, int &m, int &n, int& kL, int& kU, float& alpha, float *A, int &ldA, float *x, int& incx, float &beta, float *y, int& incy)
{
   int info=-gbmv<float>(trans,m,n,kL,kU,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("SGBMV ",info);
   return 0;
}
