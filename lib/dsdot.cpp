//
//  dsdot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dotx.h"

using LATL::DOTX;

double dsdot_(int &n, float *x, int& incx, float *y, int& incy)
{
   return DOTX<float,double>(n,x,incx,y,incy);
}
