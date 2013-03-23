//
//  sspr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "spr.h"

using LATL::spr;

int sspr2_(char& uplo, int &n, float& alpha, float *x, int& incx, float *y, int& incy, float *A)
{
   int info=-spr<float>(uplo,n,alpha,x,incx,y,incy,A);
   if(info>0)
      xerbla_("SSPR2 ",info);
   return 0;
}
