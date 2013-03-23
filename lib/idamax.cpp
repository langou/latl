//
//  idamax.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "imax.h"

using LATL::imax;

int idamax_(int &n, double *x, int& incx)
{
   return imax<double>(n,x,incx)+1;
}
