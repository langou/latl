//
//  dasum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "asum.h"
using latl::asum;

double dasum_(int &n, double *x, int& incx)
{
   return asum<double>(n,x,incx);
}

