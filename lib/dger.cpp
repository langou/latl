//
//  dger.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "ger.h"

using latl::ger;

int dger_(int &m, int &n, double &alpha, double *x, int& incx, double *y, int& incy, double *A, int &ldA)
{
   int info=-ger<double>(m,n,alpha,x,incx,y,incy,A,ldA);
   if(info!=0)
      xerbla_("DGER  ",info);
   return 0;
}
