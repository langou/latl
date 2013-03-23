//
//  sdsdot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "dot.h"

using LATL::DOTX;

float sdsdot_(int &n, float *B,float *x, int& incx, float *y, int& incy)
{
   return DOTX<float,double>(n,*B,x,incx,y,incy);
}

