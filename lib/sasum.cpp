//
//  sasum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "asum.h"
using latl::asum;

float sasum_(int &n, float *x, int& incx)
{
   return asum<float>(n,x,incx);
}

