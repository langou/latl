//
//  zaxpy.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "axpy.h"

using LATL::AXPY;

int zaxpy_(int &n, complex<double>& alpha, complex<double> *x, int& incx, complex<double> *y, int& incy)
{
   AXPY<double>(n,alpha,x,incx,y,incy);
   return 0;
}
