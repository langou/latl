//
//  ddot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using LATL::dot;

double ddot_(int &n, double *x, int& incx, double *y, int& incy)
{
   return dot<double>(n,x,incx,y,incy);
}

