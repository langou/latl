//
//  ssyr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr.h"

using LATL::syr;

int ssyr_(char& uplo, int &n, float& alpha, float *x, int& incx, float *A, int &ldA)
{
   int info=-syr<float>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("SSYR  ",info);
   return 0;
}
