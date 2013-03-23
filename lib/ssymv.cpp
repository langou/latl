//
//  ssymv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "symv.h"

using LATL::SYMV;

int ssymv_(char& uplo, int &n, float& alpha, float *A, int &ldA, float *x, int& incx, float &beta, float *y, int& incy)
{
   int info=-SYMV<float>(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info!=0)
      xerbla_("SSYMV ",info);
   return 0;
}
