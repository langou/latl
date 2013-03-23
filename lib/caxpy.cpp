//
//  caxpy.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "axpy.h"

using LATL::AXPY;

int caxpy_(int &n, complex<float>& alpha, complex<float> *x, int& incx, complex<float> *y, int& incy)
{
   AXPY<float>(n,alpha,x,incx,y,incy);
   return 0;
}
