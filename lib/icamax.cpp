//
//  icamax.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "imax.h"

using LATL::IMAX;

int icamax_(int &n, complex<float> *x, int& incx)
{
   return IMAX<float>(n,x,incx)+1;
}
