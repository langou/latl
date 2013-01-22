//
//  sspmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "spmv.h"

using latl::spmv;

int sspmv_(char& uplo, int &n, float& alpha, float *A, float *x, int& incx, float &beta, float *y, int& incy)
{
   int info=-spmv<float>(uplo,n,alpha,A,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("SSPMV ",info);
   return 0;
}
