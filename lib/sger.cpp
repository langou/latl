//
//  sger.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "ger.h"

using latl::ger;

int sger_(int &m, int &n, float& alpha, float *x, int& incx, float *y, int& incy, float *A, int &ldA)
{
   int info=-ger<float>(m,n,alpha,x,incx,y,incy,A,ldA);
   if(info!=0)
      xerbla_("SGER  ",info);
   return 0;
}
