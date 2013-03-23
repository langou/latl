//
//  cdotc.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using LATL::DOTC;

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
complex<float> cdotc_(int &n, complex<float> *x, int& incx, complex<float> *y, int& incy)
{
   return DOTC<float>(n,x,incx,y,incy);
}
#else
int cdotc_(complex<float> &ret, int &n, complex<float> *x, int &incx, complex<float> *y, int &incy)
{
   ret=DOTC<float>(n,x,incx,y,incy);
   return 0;
}
#endif
