//
//  ssyr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr2.h"

using LATL::SYR2;

int ssyr2_(char& uplo, int &n, float& alpha, float *x, int& incx, float *y, int& incy, float *A, int &ldA)
{
   int info=-SYR2<float>(uplo,n,alpha,x,incx,y,incy,A,ldA);
   if(info>0)
      xerbla_("SSYR2 ",info);
   return 0;
}
