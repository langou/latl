//
//  csrot.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rot.h"

using LATL::ROT;

int csrot_(int &n, complex<float> *x, int& incx, complex<float> *y, int& incy, float& C, float& S)
{
   ROT<float>(n,x,incx,y,incy,C,S);
   return 0;
}

