//
//  daxpy.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "axpy.h"

using latl::axpy;

int daxpy_(int &n, double& alpha, double *x, int& incx, double *y, int& incy)
{
   axpy<double>(n,alpha,x,incx,y,incy);
   return 0;
}
