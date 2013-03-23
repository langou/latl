//
//  ssbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "sbmv.h"

using LATL::SBMV;

int ssbmv_(char& uplo, int &n, int& k, float& alpha, float *A, int &ldA, float *x, int& incx, float &beta, float *y, int& incy)
{
   int info=-SBMV<float>(uplo,n,k,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("SSBMV ",info);
   return 0;
}
