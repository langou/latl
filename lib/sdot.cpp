//
//  sdot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using LATL::dot;

float sdot_(int &n, float *x, int& incx, float *y, int& incy)
{
   return dot<float>(n,x,incx,y,incy);
}

