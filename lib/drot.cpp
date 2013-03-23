//
//  drot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rot.h"

using LATL::rot;

int drot_(int &n, double *x, int& incx, double *y, int& incy, double& C, double& S)
{
   rot<double>(n,x,incx,y,incy,C,S);
   return 0;
}

