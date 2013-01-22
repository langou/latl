//
//  zdrot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rot.h"

using latl::rot;

int zdrot_(int &n, complex<double> *x, int& incx, complex<double> *y, int& incy, double &C, double &S)
{
   rot<double>(n,x,incx,y,incy,C,S);
   return 0;
}

