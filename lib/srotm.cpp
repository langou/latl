//
//  srotm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotm.h"

using LATL::ROTM;

int srotm_(int&n , float *x, int& incx, float *y, int& incy, float *param)
{
   float flag=param[0];
   ROTM<float>(n,x,incx,y,incy,flag,param+1);
   return 0;
}
