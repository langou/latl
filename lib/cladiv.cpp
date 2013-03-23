//
//  cladiv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "ladiv.h"

using LATL::ladiv;

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
complex<float> zladiv_(complex<float> &x, complex<float> &y)
{
   return ladiv<float>(x,y);
}
#else
int cladiv_(complex<float>& ret, complex<float> &x, complex<float> &y)
{
   ret=ladiv<float>(x,y);
   return 0;
}
#endif


