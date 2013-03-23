//
//  csscal.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "scal.h"

using LATL::SCAL;

int csscal_(int &n, float& alpha, complex<float> *x, int& incx)
{
   SCAL<float>(n,alpha,x,incx);
   return 0;
}

