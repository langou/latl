//
//  dznrm2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "nrm2.h"

using LATL::nrm2;

double dznrm2_(int &n, complex<double> *x, int& incx)
{
   return nrm2<double>(n,x,incx);
}

