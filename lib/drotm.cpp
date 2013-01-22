//
//  drotm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotm.h"

using latl::rotm;

int drotm_(int&n , double *x, int& incx, double *y, int& incy, double *param)
{
   double flag=param[0];
   rotm<double>(n,x,incx,y,incy,flag,param+1);
   return 0;
}
