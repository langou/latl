//
//  zcopy.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "copy.h"

using latl::copy;

int zcopy_(int &n, complex<double> *x, int& incx, complex<double> *y, int& incy)
{
   copy<double>(n,x,incx,y,incy);
   return 0;
}

