//
//  scasum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "asum.h"

using LATL::ASUM;

float scasum_(int &n, complex<float> *x, int& incx)
{
   return ASUM<float>(n,x,incx);
}
