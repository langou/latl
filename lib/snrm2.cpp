//
//  snrm2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "nrm2.h"

using latl::nrm2;

float snrm2_(int &n, float *x, int& incx)
{
   return nrm2<float>(n,x,incx);
}
