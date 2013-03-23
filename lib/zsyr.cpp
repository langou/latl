//
//  zsyr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr.h"

using LATL::syr;

int zsyr_(char& uplo, int &n, complex<float>& alpha, complex<float> *x, int& incx, complex<float> *A, int &ldA)
{
   int info=-syr<float>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("ZSYR  ",info);
   return 0;
}
