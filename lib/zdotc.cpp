//
//  zdotc.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using latl::dotc;

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
complex<double> zdotc_(int &n, complex<double> *x, int& incx, complex<double> *y, int& incy)
{
   return dotc<double>(n,x,incx,y,incy);
}
#else
int zdotc_(complex<double> &ret, int &n, complex<double> *x, int &incx, complex<double> *y, int &incy)
{
   ret=dotc<double>(n,x,incx,y,incy);
   return 0;
}
#endif


