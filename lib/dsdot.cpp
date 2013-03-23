//
//  dsdot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using LATL::dotx;

double dsdot_(int &n, float *x, int& incx, float *y, int& incy)
{
   return dotx<float,double>(n,x,incx,y,incy);
}
